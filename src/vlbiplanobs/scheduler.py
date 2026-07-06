# -*- coding: utf-8 -*-
# Licensed under GPLv3+ - see LICENSE
"""VLBI Observation Scheduler.

Arranges scan blocks across a VLBI observation: fringe finders, polarisation
calibrators, eMERLIN 3C286, and science targets. Generates SCHED .key files
with ``group N rep R`` syntax and Jb1 source-change mitigation.

Classes
-------
ScheduledScanBlock
    Dataclass for a scan block placed at a specific time.
ObservationScheduler
    Main scheduler that orchestrates placement and key-file generation.
"""

from dataclasses import dataclass, field
from typing import Optional
from importlib import resources
from datetime import datetime
import logging
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from .sources import Source, SourceType, ScanBlock, Scan
from .observation import Observation
from . import calibrators

try:
    from ortools.sat.python import cp_model
    _HAS_ORTOOLS = True
except ImportError:  # pragma: no cover - exercised only when ortools missing
    cp_model = None  # type: ignore
    _HAS_ORTOOLS = False

log = logging.getLogger(__name__)

_EMERLIN_CODES = {'CM', 'KN', 'PI', 'DA', 'DE', 'JB2', 'JB1'}
_POLCAL_NAMES = ['3C84', 'OQ208', 'DA193']
_3C286_COORD = '13h31m08.288s +30d30m32.96s'
_JB1_MAX_SRC_CHANGES_PER_HOUR = 12


def _fmt_dur(q: u.Quantity) -> str:
    """Format an astropy duration quantity as 'M:SS' for SCHED key files.

    Parameters
    ----------
    q : Quantity
        Duration quantity.

    Returns
    -------
    str
        Formatted duration string.
    """
    total_sec = int(round(q.to(u.s).value))
    return f"{total_sec // 60}:{total_sec % 60:02d}"


def _intent_str(stype: SourceType) -> str:
    """Map a SourceType to a SCHED intent string (empty if none applies).

    Parameters
    ----------
    stype : SourceType
        Source type to map.

    Returns
    -------
    str
        SCHED intent string.
    """
    mapping = {
        SourceType.FRINGEFINDER: 'FRINGE_FINDER',
        SourceType.AMPLITUDECAL: 'EMERLIN_AMP',
        SourceType.POLCAL: 'POLCAL',
        SourceType.CHECKSOURCE: 'CHECK',
        SourceType.TARGET: 'TARGET',
        SourceType.PHASECAL: 'PHASE_CAL',
    }
    return mapping.get(stype, '')


# ======================================================================
# Data classes
# ======================================================================

@dataclass
class ScheduledScanBlock:
    """A scan block scheduled at a specific time with timing and quality metrics.

    Attributes
    ----------
    name : str
        Identifier (e.g. 'FF_1', 'R1_D').
    block : ScanBlock
        The original ScanBlock being scheduled.
    start_time : Time
        Start time of this scheduled block.
    end_time : Time
        End time of this scheduled block.
    scans : list[Scan]
        Expanded list of scans filling this time slot.
    n_antennas : int
        Number of antennas that can observe during this block.
    mean_elevation : float
        Mean elevation across observing antennas (degrees).
    """
    name: str
    block: ScanBlock
    start_time: Time
    end_time: Time
    scans: list[Scan] = field(default_factory=list)
    n_antennas: int = 0
    mean_elevation: float = 0.0

    @property
    def duration(self) -> u.Quantity:
        """Duration of this scheduled block."""
        return (self.end_time - self.start_time).to(u.min)


# ======================================================================
# Scheduler
# ======================================================================

class ObservationScheduler:
    """Schedule VLBI scan blocks across the observation time range.

    Parameters
    ----------
    observation : Observation
        The observation containing scan blocks to schedule.
    min_antennas : int
        Minimum antennas required for a valid time slot.
    require_all_antennas : bool
        If True, only schedule when all antennas can observe.
    fringefinder_spec : list[str] or None
        Source names, or ``['N']`` to auto-select *N* FF sources.
    polcal : bool
        Whether polarisation calibration scans are required.
    """

    FF_DUR = 5 * u.min
    FF_INTERVAL = 2 * u.h
    POLCAL_DUR = 5 * u.min
    EMERLIN_3C286_DUR = 5 * u.min
    GAP_DUR = 30 * u.s
    GAP_INTERVAL = 10 * u.min

    # --- CP-SAT objective weights (integer coefficients) ---
    # Spread term rewards a large minimum gap (in grid cells) between
    # consecutive fringe-finder scans.  Separation term penalises the angular
    # distance (deg) between a fringe finder and the science target it
    # interrupts.  Elevation term mildly prefers high-elevation FF placements.
    CP_W_SPREAD = 1000
    CP_W_SEPARATION = 10
    CP_W_ELEVATION = 1
    CP_MAX_SOLVE_SECONDS = 5.0

    def __init__(self, observation: Observation, min_antennas: int = 2,
                 require_all_antennas: bool = False,
                 fringefinder_spec: Optional[list[str]] = None, polcal: bool = False):
        self.obs = observation
        self.min_ant = min_antennas
        self.require_all = require_all_antennas
        self._ff_spec = fringefinder_spec or ['2']
        self._ff_n_scans: Optional[int] = None
        if len(self._ff_spec) == 1 and self._ff_spec[0].isdigit():
            self._ff_n_scans = int(self._ff_spec[0])
        self._polcal = polcal
        self._scheduled: list[ScheduledScanBlock] = []
        self._precompute()

    # ------------------------------------------------------------------
    # Pre-computation
    # ------------------------------------------------------------------

    def _precompute(self):
        """Build per-block visibility and mean-elevation arrays over obs.times.

        Precomputes visibility counts and mean elevations for all scan blocks
        at all observation time steps for efficient scheduling.
        """
        self._blocks = list(self.obs.scans.keys())
        self._n_times = len(self.obs.times)
        self._n_ant = len(self.obs.stations)
        self._dt = (self.obs.times[1] - self.obs.times[0]).to(u.min).value
        is_obs = self.obs.is_observable()
        elevs = self.obs.elevations()
        self._vis: dict[str, np.ndarray] = {}
        self._elev: dict[str, np.ndarray] = {}
        self._is_ff: dict[str, bool] = {}
        for name in self._blocks:
            block = self.obs.scans[name]
            self._is_ff[name] = block.has(SourceType.FRINGEFINDER)
            self._vis[name] = (np.sum(np.array(list(is_obs[name].values())), axis=0)
                               if name in is_obs else np.zeros(self._n_times, dtype=int))
            targets = block.sources(SourceType.TARGET) or block.sources(SourceType.FRINGEFINDER) or block.sources()
            if targets and targets[0].name in elevs:
                self._elev[name] = np.nanmean(
                    [elevs[targets[0].name][a].value for a in elevs[targets[0].name]], axis=0)
            else:
                self._elev[name] = np.zeros(self._n_times)

    def _ff_blocks(self) -> list[str]:
        """Get list of fringe-finder block names.

        Returns
        -------
        list[str]
            Names of blocks containing fringe-finder sources.
        """
        return [n for n in self._blocks if self._is_ff.get(n, False)]

    def _sci_blocks(self) -> list[str]:
        """Get list of science block names.

        Returns
        -------
        list[str]
            Names of blocks not containing fringe-finder sources.
        """
        return [n for n in self._blocks if not self._is_ff.get(n, False)]

    def _t2i(self, t: Time) -> int:
        """Convert a Time to the nearest index in obs.times.

        Parameters
        ----------
        t : Time
        Time to convert.

        Returns
        -------
        int
        Index in obs.times array.
        """
        return int(np.argmin(np.abs(self.obs.times.mjd - t.mjd)))

    def _i2t(self, i: int) -> Time:
        """Convert an index to a Time from obs.times.

        Parameters
        ----------
        i : int
        Index in obs.times array.

        Returns
        -------
        Time
        Time at the given index.
        """
        return self.obs.times[min(i, self._n_times - 1)]

    def _vis_at(self, name: str, t: Time) -> int:
        """Get number of visible antennas for a block at a specific time.

        Parameters
        ----------
        name : str
        Block name.
        t : Time
        Observation time.

        Returns
        -------
        int
        Number of visible antennas.
        """
        return int(self._vis[name][self._t2i(t)])

    def _elev_at(self, name: str, t: Time) -> float:
        """Get mean elevation for a block at a specific time.

        Parameters
        ----------
        name : str
        Block name.
        t : Time
        Observation time.

        Returns
        -------
        float
        Mean elevation in degrees.
        """
        return float(self._elev[name][self._t2i(t)])

    def _vis_mean(self, name: str, t0: Time, t1: Time) -> float:
        """Mean number of antennas visible for a block between two times.

        Parameters
        ----------
        name : str
        Block name.
        t0 : Time
        Start time.
        t1 : Time
        End time.

        Returns
        -------
        float
        Mean number of visible antennas.
        """
        i0, i1 = self._t2i(t0), self._t2i(t1)
        arr = self._vis[name][i0:max(i1, i0 + 1)]
        return float(np.mean(arr))

    def _elev_mean(self, name: str, t0: Time, t1: Time) -> float:
        """Mean elevation for a block between two times.

        Parameters
        ----------
        name : str
        Block name.
        t0 : Time
        Start time.
        t1 : Time
        End time.

        Returns
        -------
        float
        Mean elevation in degrees.
        """
        i0, i1 = self._t2i(t0), self._t2i(t1)
        arr = self._elev[name][i0:max(i1, i0 + 1)]
        return float(np.nanmean(arr))

    # ------------------------------------------------------------------
    # Slot helpers
    # ------------------------------------------------------------------

    def _find_best(self, name: str, t0: Time, t1: Time,
                   dur: u.Quantity) -> Optional[tuple[Time, Time, int, float]]:
        """Find optimal placement for a block in a time range with given duration.

        Parameters
        ----------
        name : str
        Block name.
        t0 : Time
        Start of search window.
        t1 : Time
        End of search window.
        dur : Quantity
        Required duration.

        Returns
        -------
        tuple[Time, Time, int, float] or None
        (start, end, min_antennas, mean_elevation) or None if no valid placement.
        """
        i0, i1 = self._t2i(t0), self._t2i(t1)
        steps = max(1, int(np.ceil(dur.to(u.min).value / self._dt)))
        if i1 - i0 < steps:
            return None
        vis, elev = self._vis[name], self._elev[name]
        min_req = self._n_ant if self.require_all else self.min_ant
        best_score, best_i = -1.0, -1
        for i in range(i0, i1 - steps + 1):
            n = int(np.min(vis[i:i + steps]))
            e = float(np.mean(elev[i:i + steps]))
            if n >= min_req and not np.isnan(e):
                score = n * 100.0 + e
                if score > best_score:
                    best_score, best_i = score, i
        if best_i < 0:
            return None
        return (self._i2t(best_i), self._i2t(best_i + steps),
                int(np.min(vis[best_i:best_i + steps])),
                float(np.mean(elev[best_i:best_i + steps])))

    def _get_slots(self, sched: list[ScheduledScanBlock], t0: Time, t1: Time) -> list[tuple[Time, Time]]:
        """Return free time slots between already-scheduled blocks.

        Parameters
        ----------
        sched : list[ScheduledScanBlock]
        Already scheduled blocks.
        t0 : Time
        Observation start time.
        t1 : Time
        Observation end time.

        Returns
        -------
        list[tuple[Time, Time]]
        List of (start, end) time tuples for free slots.
        """
        if not sched:
            return [(t0, t1)]
        s = sorted(sched, key=lambda b: b.start_time.mjd)
        slots: list[tuple[Time, Time]] = []
        if s[0].start_time > t0:
            slots.append((t0, s[0].start_time))
        for i in range(len(s) - 1):
            if s[i + 1].start_time > s[i].end_time:
                slots.append((s[i].end_time, s[i + 1].start_time))
        if s[-1].end_time < t1:
            slots.append((s[-1].end_time, t1))
        return slots

    def _register_block(self, name: str):
        """Register a dynamically-added block in pre-computed arrays.

        Parameters
        ----------
        name : str
        Block name to register.
        """
        if name in self._vis:
            return
        self._blocks.append(name)
        block = self.obs.scans[name]
        self._is_ff[name] = block.has(SourceType.FRINGEFINDER)
        self._vis[name] = np.full(self._n_times, self._n_ant, dtype=int)
        self._elev[name] = np.full(self._n_times, 45.0)

    # ------------------------------------------------------------------
    # Fringe-finder selection  (NEW)
    # ------------------------------------------------------------------

    def _select_fringefinders(self) -> list[str]:
        """Select FF sources and register them as scan blocks.

        Strategy
        --------
        1. If the user gave explicit source names, use those.
        2. Otherwise find FFs visible by ALL antennas for the ENTIRE observation;
           pick the one closest to the mean science-target position.
        3. Fallback: find FFs visible by at least *min_ant* antennas; for each
           FF time slot pick the best-visible source closest to the targets.

        Returns
        -------
        list[str]
            Block names registered in ``self.obs.scans``.
        """
        existing = self._ff_blocks()
        if existing:
            return existing

        spec = self._ff_spec
        if spec and not (len(spec) == 1 and spec[0].isdigit()):
            return self._register_ff_sources(self._lookup_ff_sources(spec))

        # Try all-antenna, all-time visibility first
        cands, _, _ = calibrators.get_fringe_finder_sources(
            self.obs.stations, self.obs.times,
            min_elevation=20 * u.deg, min_flux=0.5 * u.Jy, require_all_stations=True)

        if cands:
            best = self._pick_closest_ff(cands)
            return self._register_ff_sources([best])

        # Relax to partial visibility
        cands, _, _ = calibrators.get_fringe_finder_sources(
            self.obs.stations, self.obs.times,
            min_elevation=20 * u.deg, min_flux=0.5 * u.Jy, require_all_stations=False)

        if cands:
            best = self._pick_closest_ff(cands)
            return self._register_ff_sources([best])

        log.warning("No fringe-finder sources found for this observation.")
        return []

    def _pick_closest_ff(self, candidates: list) -> Source:
        """From a list of CalibratorSource candidates, pick the one closest to science targets.

        Parameters
        ----------
        candidates : list
        List of CalibratorSource candidates.

        Returns
        -------
        Source
        The closest source converted to a Source object.
        """
        mean_coord = self._mean_target_coord()
        if mean_coord is None:
            src = candidates[0]
            return Source(name=src.name, coordinates=src.coord,
                          source_type=SourceType.FRINGEFINDER, other_names=[src.ivsname])
        best_src, best_sep = candidates[0], float('inf')
        for c in candidates:
            sep = c.coord.separation(mean_coord).deg
            if sep < best_sep:
                best_sep, best_src = sep, c
        return Source(name=best_src.name, coordinates=best_src.coord,
                      source_type=SourceType.FRINGEFINDER, other_names=[best_src.ivsname])

    def _mean_target_coord(self) -> Optional[SkyCoord]:
        """Return the mean sky position of all science targets, or None.

        Returns
        -------
        SkyCoord or None
        Mean coordinate of all science targets, or None if no targets.
        """
        ras, decs = [], []
        for name in self._sci_blocks():
            for src in self.obs.scans[name].sources(SourceType.TARGET):
                ras.append(src.coord.ra.deg)
                decs.append(src.coord.dec.deg)
        if not ras:
            return None
        return SkyCoord(float(np.mean(ras)), float(np.mean(decs)), unit=u.deg)

    def _register_ff_sources(self, sources: list[Source]) -> list[str]:
        """Register unique FF Source objects as scan blocks.

        Parameters
        ----------
        sources : list[Source]
        List of Source objects to register.

        Returns
        -------
        list[str]
        Block names registered in obs.scans.
        """
        seen: set[str] = set()
        names: list[str] = []
        for src in sources:
            if src.name in seen:
                continue
            seen.add(src.name)
            block_name = f"FF_{src.name}"
            if block_name not in self.obs.scans:
                self.obs.scans[block_name] = ScanBlock([Scan(src, duration=self.FF_DUR)])
                self._register_block(block_name)
            names.append(block_name)
        return names

    def _lookup_ff_sources(self, names: list[str]) -> list[Source]:
        """Look up named fringe-finder sources from the RFC catalog.

        Each entry in *names* can be 'name/coordinates' to supply explicit
        coordinates, plain coordinates, or a catalog name.

        Parameters
        ----------
        names : list[str]
        List of source specifications.

        Returns
        -------
        list[Source]
        List of Source objects with FRINGEFINDER type.
        """
        cat = calibrators.RFCCatalog(min_flux=0.0 * u.Jy, band='c')
        result: list[Source] = []
        for spec in names:
            parsed_name, parsed_coord = Source.parse_source_spec(spec)
            if parsed_name is not None and parsed_coord is not None:
                result.append(Source(
                    parsed_name, coordinates=Source._parse_coord_str(parsed_coord),
                    source_type=SourceType.FRINGEFINDER))
            else:
                lookup = parsed_name or spec
                rfc_src = cat.get_source(lookup)
                if rfc_src is not None:
                    result.append(Source(name=rfc_src.name, coordinates=rfc_src.coord,
                                        source_type=SourceType.FRINGEFINDER,
                                        other_names=[rfc_src.ivsname]))
                else:
                    try:
                        result.append(Source.source_from_str(spec, source_type=SourceType.FRINGEFINDER))
                    except ValueError:
                        log.warning("Fringe-finder source '%s' not found — skipping.", spec)
        return result

    # ------------------------------------------------------------------
    # Polcal / eMERLIN helpers
    # ------------------------------------------------------------------

    def _create_polcal_blocks(self) -> list[str]:
        """Create 3C84 / OQ208 / DA193 polcal blocks.

        Returns
        -------
        list[str]
            Block names registered in obs.scans.
        """
        cat = calibrators.RFCCatalog(min_flux=0.0 * u.Jy, band='c')
        added: list[str] = []
        for polcal_name in _POLCAL_NAMES:
            rfc_src = cat.get_source(polcal_name)
            if rfc_src is None:
                continue
            src = Source(name=rfc_src.name, coordinates=rfc_src.coord,
                         source_type=SourceType.POLCAL, other_names=[rfc_src.ivsname])
            block_name = f"POLCAL_{rfc_src.ivsname}"
            self.obs.scans[block_name] = ScanBlock([Scan(src, duration=self.POLCAL_DUR)])
            self._register_block(block_name)
            added.append(block_name)
        return added

    def _has_emerlin(self) -> bool:
        """Check if eMERLIN stations are in the array.

        Returns
        -------
        bool
        True if any eMERLIN station is present.
        """
        return any(s.codename.upper() in _EMERLIN_CODES for s in self.obs.stations)

    def _has_jb1(self) -> bool:
        """Return True if Jodrell Bank Mk2 (Jb1 / JB1) is in the array.

        Returns
        -------
        bool
        True if JB1 is present.
        """
        return any(s.codename.upper() == 'JB1' for s in self.obs.stations)

    def _create_emerlin_3c286_block(self) -> Optional[str]:
        """Create a 3C286 scan block for eMERLIN flux-scale calibration.

        Returns
        -------
        str or None
        Block name if created, None otherwise.
        """
        src = Source(name='3C286', coordinates=_3C286_COORD, source_type=SourceType.AMPLITUDECAL)
        block_name = 'eMERLIN_3C286'
        self.obs.scans[block_name] = ScanBlock([Scan(src, duration=self.EMERLIN_3C286_DUR)])
        self._register_block(block_name)
        return block_name

    # ------------------------------------------------------------------
    # FF scheduling  (REWRITTEN — evenly spread, never consecutive)
    # ------------------------------------------------------------------

    def _effective_obs_window(self, t0: Time, t1: Time) -> tuple[Time, Time]:
        """Return the effective observation window where science is possible.

        Shrinks [t0, t1] to the range where at least one science block has
        >= min_ant antennas visible.  FF scans should be placed only within
        this window so they bracket actual science time.

        Parameters
        ----------
        t0 : Time
        Observation start time.
        t1 : Time
        Observation end time.

        Returns
        -------
        tuple[Time, Time]
        (effective_start, effective_end).
        """
        sci_names = self._sci_blocks()
        if not sci_names:
            return t0, t1
        # Combine visibility across all science blocks (take max per time step)
        combined = np.zeros(self._n_times, dtype=int)
        for n in sci_names:
            combined = np.maximum(combined, self._vis[n])
        ok = np.where(combined >= self.min_ant)[0]
        if len(ok) == 0:
            return t0, t1
        eff_t0 = self._i2t(int(ok[0]))
        eff_t1 = self._i2t(int(ok[-1]))
        # Pad the end by FF_DUR so we can fit a closing FF scan
        eff_t1 = min(t1, eff_t1 + self.FF_DUR)
        return max(t0, eff_t0), eff_t1

    def _schedule_ff(self, ff_names: list[str], t0: Time, t1: Time) -> list[ScheduledScanBlock]:
        """Place FF scans evenly across the observation, never adjacent.

        One scan at start, one at end, additional every ~FF_INTERVAL in between.
        FF scans are only placed within the effective observation window where
        science targets are observable.

        Parameters
        ----------
        ff_names : list[str]
            Block names for FF sources (typically one element).
        t0 : Time
            Observation start time.
        t1 : Time
            Observation end time.

        Returns
        -------
        list[ScheduledScanBlock]
            Scheduled FF blocks in chronological order.
        """
        if not ff_names:
            return []

        # Restrict to where science is actually observable
        eff_t0, eff_t1 = self._effective_obs_window(t0, t1)
        dur = self.FF_DUR
        obs_dur_h = (eff_t1 - eff_t0).to(u.h).value

        if self._ff_n_scans is not None:
            n_ff = max(1, self._ff_n_scans)
            log.info("Using user-requested FF scan count: %d", n_ff)
        elif obs_dur_h <= 1.5:
            n_ff = 1
        elif obs_dur_h <= 3.0:
            n_ff = 2
        else:
            n_mid = max(0, int((obs_dur_h - dur.to(u.h).value) / self.FF_INTERVAL.to(u.h).value))
            n_ff = n_mid + 1

        total_ff_time = n_ff * dur
        total_sci_time = (eff_t1 - eff_t0) - total_ff_time
        if total_sci_time.to(u.min).value < 0:
            total_sci_time = 0 * u.min

        gap = total_sci_time / max(n_ff - 1, 1) if n_ff > 1 else 0 * u.min
        starts: list[Time] = []
        t = eff_t0
        for _ in range(n_ff):
            starts.append(t)
            t = t + dur + gap

        result: list[ScheduledScanBlock] = []
        for idx, ts in enumerate(starts):
            te = ts + dur
            block_name = ff_names[idx % len(ff_names)]
            n_vis = self._vis_at(block_name, ts)
            if n_vis < self.min_ant:
                log.info("Skipping FF slot at %s — only %d antennas can observe.", ts.iso, n_vis)
                continue
            block = self.obs.scans[block_name]
            result.append(ScheduledScanBlock(
                name=f"FF_{len(result) + 1}", block=block, start_time=ts, end_time=te,
                scans=list(block.scans), n_antennas=n_vis,
                mean_elevation=self._elev_at(block_name, ts)))
        return result

    # ------------------------------------------------------------------
    # Polcal / single-block scheduling  (unchanged logic)
    # ------------------------------------------------------------------

    def _schedule_polcal(self, polcal_names: list[str],
                         slots: list[tuple[Time, Time]]) -> list[ScheduledScanBlock]:
        """Spread polcal scans at 10 %, 50 %, 90 % of available time.

        Parameters
        ----------
        polcal_names : list[str]
            Block names for polcal sources.
        slots : list[tuple[Time, Time]]
            Available time slots.

        Returns
        -------
        list[ScheduledScanBlock]
            Scheduled polcal blocks.
        """
        if not polcal_names or not slots:
            return []
        n_polcal = min(3, len(polcal_names))
        dur = self.POLCAL_DUR
        total = sum(((s1 - s0).to(u.min).value for s0, s1 in slots)) * u.min
        if total < dur:
            return []
        placed: list[ScheduledScanBlock] = []
        for frac, pc_name in zip([0.1, 0.5, 0.9][:n_polcal], polcal_names):
            block = self.obs.scans[pc_name]
            elapsed = 0.0 * u.min
            target_time = frac * total
            for s0, s1 in slots:
                slot_dur = (s1 - s0).to(u.min)
                if elapsed + slot_dur >= target_time and slot_dur >= dur:
                    offset = max(target_time - elapsed, 0.0 * u.min)
                    ts = s0 + offset
                    if ts + dur <= s1:
                        placed.append(ScheduledScanBlock(
                            name=f"POLCAL_{pc_name.split('_')[-1]}", block=block,
                            start_time=ts, end_time=ts + dur, scans=[block.scans[0]],
                            n_antennas=self._vis_at(pc_name, ts),
                            mean_elevation=self._elev_at(pc_name, ts)))
                        break
                elapsed += slot_dur
        return placed

    def _schedule_single_block(self, block_name: str, slots: list[tuple[Time, Time]],
                               label: str) -> Optional[ScheduledScanBlock]:
        """Place a single auxiliary block in the best available slot.

        Parameters
        ----------
        block_name : str
        Block name to schedule.
        slots : list[tuple[Time, Time]]
        Available time slots.
        label : str
        Label for the scheduled block.

        Returns
        -------
        ScheduledScanBlock or None
        Scheduled block or None if no valid slot found.
        """
        block = self.obs.scans[block_name]
        dur = sum(s.duration.to(u.min).value for s in block.scans) * u.min
        best = None
        for s0, s1 in slots:
            if (s1 - s0).to(u.min) >= dur:
                r = self._find_best(block_name, s0, s1, dur)
                if r and (not best or r[2] * 100 + r[3] > best[2] * 100 + best[3]):
                    best = r
        if best:
            ts, te, n, e = best
            return ScheduledScanBlock(name=label, block=block, start_time=ts, end_time=te,
                                      scans=[block.scans[0]], n_antennas=n, mean_elevation=e)
        return None

    # ------------------------------------------------------------------
    # Science scheduling  (REWRITTEN — optimise per slot)
    # ------------------------------------------------------------------

    def _schedule_sci(self, names: list[str], slots: list[tuple[Time, Time]]) -> list[ScheduledScanBlock]:
        """Fill each free slot with the science block that has the best observing conditions there.

        For a single science block all slots are filled with it.  For multiple
        blocks the one with the highest (n_antennas * 100 + mean_elevation)
        score wins each slot.

        Parameters
        ----------
        names : list[str]
            Science block names.
        slots : list[tuple[Time, Time]]
            Available time windows between FF / calibrator scans.

        Returns
        -------
        list[ScheduledScanBlock]
            Scheduled science blocks.
        """
        if not names or not slots:
            return []
        scheduled: list[ScheduledScanBlock] = []
        for s0, s1 in slots:
            slot_dur = (s1 - s0).to(u.min)
            best_name: Optional[str] = None
            best_score = -1.0
            for name in names:
                block = self.obs.scans[name]
                min_dur = sum(s.duration.to(u.min).value for s in block.scans) * u.min
                if slot_dur < min_dur:
                    continue
                score = self._vis_mean(name, s0, s1) * 100.0 + self._elev_mean(name, s0, s1)
                if score > best_score:
                    best_score, best_name = score, name
            if best_name is None:
                continue
            block = self.obs.scans[best_name]
            try:
                scans = block.fill(slot_dur)
            except ValueError:
                scans = list(block.scans)
            mid = s0 + (s1 - s0) * 0.5
            scheduled.append(ScheduledScanBlock(
                name=best_name, block=block, start_time=s0, end_time=s1, scans=scans,
                n_antennas=self._vis_at(best_name, mid),
                mean_elevation=self._elev_at(best_name, mid)))
        return scheduled

    # ------------------------------------------------------------------
    # Multi-source scheduling helpers
    # ------------------------------------------------------------------

    def _source_main_coord(self, name: str) -> Optional[SkyCoord]:
        """Return the sky coordinate of the primary science target in a scan block.

        Parameters
        ----------
        name : str
            Block name in obs.scans.

        Returns
        -------
        SkyCoord or None
            Coordinate of the primary target, or None if no target found.
        """
        block = self.obs.scans.get(name)
        if block is None:
            return None
        targets = block.sources(SourceType.TARGET) or block.sources(SourceType.PULSAR) or block.sources()
        return targets[0].coord if targets else None

    def _find_optimal_window(self, name: str, t0: Time, t1: Time, dur: u.Quantity) -> tuple[Time, float]:
        """Find the start time of the best observing window for a block.

        Searches the interval [t0, t1] for the placement of *dur* that maximises
        n_antennas * 100 + mean_elevation.

        Parameters
        ----------
        name : str
        Block name.
        t0 : Time
        Start of search window.
        t1 : Time
        End of search window.
        dur : Quantity
        Required duration.

        Returns
        -------
        tuple[Time, float]
        (best_start, score) — score is -1 if no valid window found.
        """
        result = self._find_best(name, t0, t1, dur)
        if result:
            ts, _te, n, e = result
            return ts, float(n) * 100.0 + e
        return t0, -1.0

    def _nearest_neighbor_order(self, names: list[str]) -> list[str]:
        """Reorder names using nearest-neighbour heuristic to minimise total angular slewing.

        Starts from the first element of *names* (which should already be sorted
        by optimal observation window) and greedily picks the geographically
        closest remaining source at each step.

        Parameters
        ----------
        names : list[str]
            Source block names, pre-sorted by optimal start time.

        Returns
        -------
        list[str]
            Reordered names.
        """
        if len(names) <= 2:
            return list(names)
        coords: dict[str, Optional[SkyCoord]] = {n: self._source_main_coord(n) for n in names}
        remaining = list(names)
        ordered = [remaining.pop(0)]
        while remaining:
            last_coord = coords.get(ordered[-1])
            if last_coord is None:
                ordered.append(remaining.pop(0))
                continue
            best_next = min(remaining, key=lambda n: (
                last_coord.separation(coords[n]).deg if coords.get(n) is not None else 360.0))
            ordered.append(best_next)
            remaining.remove(best_next)
        return ordered

    def _choose_ff_positions(self, n_sources: int, n_ff: int) -> set[int]:
        """Choose at which source boundaries to insert FF scans.

        Positions range from 0 (before source 0) to n_sources (after last source).
        FFs are spread as evenly as possible across the observation.

        Parameters
        ----------
        n_sources : int
        Number of science sources.
        n_ff : int
        Number of FF scans to place.

        Returns
        -------
        set[int]
        Set of boundary positions where FF scans should be inserted.
        """
        if n_ff <= 0:
            return set()
        if n_ff >= n_sources + 1:
            return set(range(n_sources + 1))
        positions: set[int] = set()
        for k in range(n_ff):
            idx = round((k + 0.5) * (n_sources + 1) / n_ff)
            positions.add(min(max(idx, 0), n_sources))
        return positions

    def _determine_ff_count(self, obs_dur_h: float) -> int:
        """Determine the number of FF scans for a given observation duration.

        Respects user-specified count.  Otherwise:
        <= 1.5 h → 1, <= 3 h → 2, > 3 h → 2 + one per FF_INTERVAL.

        Parameters
        ----------
        obs_dur_h : float
            Observation duration in hours.

        Returns
        -------
        int
            Number of FF scans to schedule.
        """
        if self._ff_n_scans is not None:
            return max(1, self._ff_n_scans)
        if obs_dur_h <= 1.5:
            return 1
        if obs_dur_h <= 3.0:
            return 2
        extra = max(0, int(obs_dur_h / self.FF_INTERVAL.to(u.h).value) - 1)
        return 2 + extra

    def _schedule_multi_source(self, sci_names: list[str], ff_names: list[str],
                                t0: Time, t1: Time) -> list[ScheduledScanBlock]:
        """Schedule multiple science sources with equalised time and interleaved FF scans.

        Algorithm
        ---------
        1. Calculate equal science time per source (total_time - ff_overhead).
        2. Find the optimal observing window for each source.
        3. Sort sources chronologically by their optimal window; break ties with
           nearest-neighbour sky distance to minimise slewing.
        4. Place FF scans at evenly-distributed boundaries between source blocks.
        5. Place each source block sequentially; the duration is adjusted so that
           the last block ends exactly at t1.

        Parameters
        ----------
        sci_names : list[str]
            Science block names (already filtered, no FF/polcal/eMERLIN blocks).
        ff_names : list[str]
            FF block names registered in obs.scans.
        t0 : Time
            Observation start time.
        t1 : Time
            Observation end time.

        Returns
        -------
        list[ScheduledScanBlock]
            All scheduled blocks (science + FF) in chronological order.
        """
        n_sci = len(sci_names)
        obs_dur_h = (t1 - t0).to(u.h).value
        n_ff = self._determine_ff_count(obs_dur_h)

        # --- time budget ---
        total_ff_time = n_ff * self.FF_DUR
        sci_time_total = max((t1 - t0) - total_ff_time, (t1 - t0) * 0.1)
        time_per_source = sci_time_total / n_sci

        # --- find optimal windows ---
        source_peak: dict[str, Time] = {}
        for name in sci_names:
            peak_start, _ = self._find_optimal_window(name, t0, t1, time_per_source)
            source_peak[name] = peak_start

        # --- sort by optimal window start, then minimise slewing ---
        ordered = sorted(sci_names, key=lambda n: source_peak[n].mjd)
        ordered = self._nearest_neighbor_order(ordered)

        # --- choose where to insert FF scans ---
        ff_positions = self._choose_ff_positions(n_sci, n_ff)
        log.info("Multi-source schedule: %d sources, %d FF scans at positions %s",
                 n_sci, n_ff, sorted(ff_positions))

        # --- build timeline ---
        all_blocks: list[ScheduledScanBlock] = []
        ff_count = 0
        current_t = t0

        def _place_ff(t: Time, count: int) -> Time:
            """Place one FF block starting at t; returns new current time."""
            if not ff_names:
                return t
            fn = ff_names[count % len(ff_names)]
            ff_block = self.obs.scans[fn]
            n_vis = self._vis_at(fn, t)
            all_blocks.append(ScheduledScanBlock(
                name=f"FF_{count + 1}", block=ff_block,
                start_time=t, end_time=t + self.FF_DUR,
                scans=list(ff_block.scans), n_antennas=n_vis,
                mean_elevation=self._elev_at(fn, t)))
            return t + self.FF_DUR

        for i, name in enumerate(ordered):
            # Possibly insert FF before this source
            if i in ff_positions and ff_count < n_ff:
                current_t = _place_ff(current_t, ff_count)
                ff_count += 1

            # Allocate time: divide remaining time equally among remaining sources,
            # reserving time for any FFs still to be placed after this source.
            remaining_sci_blocks = n_sci - i
            remaining_ff_time = (n_ff - ff_count) * self.FF_DUR
            available_to_end = (t1 - current_t) - remaining_ff_time
            this_dur = max(available_to_end / remaining_sci_blocks, 1.0 * u.min)

            block = self.obs.scans[name]
            sci_dur = min(this_dur, (t1 - current_t)).to(u.min)
            try:
                scans = block.fill(sci_dur)
            except (ValueError, Exception):
                scans = list(block.scans)

            # Use actual scan durations for precise time tracking (fill() may not use the full sci_dur)
            actual_dur = sum(s.duration.to(u.min).value for s in scans) * u.min
            sci_end = current_t + actual_dur
            mid = current_t + actual_dur * 0.5
            all_blocks.append(ScheduledScanBlock(
                name=name, block=block, start_time=current_t, end_time=sci_end,
                scans=scans, n_antennas=self._vis_at(name, mid),
                mean_elevation=self._elev_at(name, mid)))
            current_t = sci_end

        # Place any remaining FFs at the end
        while ff_count < n_ff and ff_names and current_t + self.FF_DUR <= t1:
            current_t = _place_ff(current_t, ff_count)
            ff_count += 1

        return sorted(all_blocks, key=lambda b: b.start_time.mjd)

    # ------------------------------------------------------------------
    # CP-SAT scheduling  (constraint-programming layout)
    # ------------------------------------------------------------------

    def _block_source_coord(self, name: str) -> Optional[SkyCoord]:
        """Return the sky coordinate of the first source in a scan block.

        Parameters
        ----------
        name : str
        Block name.

        Returns
        -------
        SkyCoord or None
        Coordinate of the first source, or None if block not found.
        """
        block = self.obs.scans.get(name)
        if block is None:
            return None
        srcs = block.sources()
        return srcs[0].coord if srcs else None

    def _min_block_duration(self, name: str) -> u.Quantity:
        """Return the minimum duration to fit one pass of a scan block.

        Parameters
        ----------
        name : str
        Block name.

        Returns
        -------
        Quantity
        Minimum duration in minutes.
        """
        block = self.obs.scans[name]
        return sum(s.duration.to(u.min).value for s in block.scans) * u.min

    def _layout_science(self, sci_names: list[str], eff_t0: Time, eff_t1: Time) -> list[ScheduledScanBlock]:
        """Lay out science blocks contiguously across [eff_t0, eff_t1] (no calibrators).

        For a single science block the whole window is filled with it.  For
        multiple blocks the time is equalised per source, ordered by optimal
        observing window, and slew-minimised via nearest-neighbour.  Fringe
        finders and other calibrators are inserted afterwards by splitting these
        science blocks (see ``_insert_aux``).

        Parameters
        ----------
        sci_names : list[str]
            Science block names.
        eff_t0 : Time
            Effective observation window start where science is observable.
        eff_t1 : Time
            Effective observation window end where science is observable.

        Returns
        -------
        list[ScheduledScanBlock]
            Contiguous science blocks covering the window.
        """
        if not sci_names:
            return []
        if len(sci_names) == 1:
            name = sci_names[0]
            block = self.obs.scans[name]
            span = (eff_t1 - eff_t0)
            try:
                scans = block.fill(span.to(u.min))
            except Exception:
                scans = list(block.scans)
            mid = eff_t0 + span * 0.5
            return [ScheduledScanBlock(
                name=name, block=block, start_time=eff_t0, end_time=eff_t1, scans=scans,
                n_antennas=self._vis_at(name, mid), mean_elevation=self._elev_at(name, mid))]

        n = len(sci_names)
        tps = (eff_t1 - eff_t0) / n
        peak = {nm: self._find_optimal_window(nm, eff_t0, eff_t1, tps)[0] for nm in sci_names}
        ordered = sorted(sci_names, key=lambda nm: peak[nm].mjd)
        ordered = self._nearest_neighbor_order(ordered)

        blocks: list[ScheduledScanBlock] = []
        cur = eff_t0
        for i, name in enumerate(ordered):
            remaining = n - i
            this_dur = (eff_t1 - cur) / remaining
            block = self.obs.scans[name]
            sci_dur = min(this_dur, (eff_t1 - cur)).to(u.min)
            try:
                scans = block.fill(sci_dur)
            except Exception:
                scans = list(block.scans)
            actual = sum(s.duration.to(u.min).value for s in scans) * u.min
            if actual.to(u.min).value <= 0:
                actual = sci_dur
            end = cur + actual
            mid = cur + actual * 0.5
            blocks.append(ScheduledScanBlock(
                name=name, block=block, start_time=cur, end_time=end, scans=scans,
                n_antennas=self._vis_at(name, mid), mean_elevation=self._elev_at(name, mid)))
            cur = end
        return blocks

    def _occupant_at(self, layout: list[ScheduledScanBlock], t: Time) -> Optional[tuple[str, Optional[SkyCoord]]]:
        """Return (block_name, target_coord) of the science block covering time t, or None.

        Parameters
        ----------
        layout : list[ScheduledScanBlock]
        Scheduled science blocks.
        t : Time
        Time to check.

        Returns
        -------
        tuple[str, SkyCoord] or None
        (block_name, target_coord) if t is within a block, else None.
        """
        for sb in layout:
            if sb.start_time <= t < sb.end_time:
                return sb.name, self._source_main_coord(sb.name)
        return None

    def _cpsat_ff_placement(self, ff_names: list[str], layout: list[ScheduledScanBlock],
                            eff_t0: Time, eff_t1: Time) -> list[tuple[Time, str]]:
        """Decide fringe-finder times and sources via CP-SAT optimisation.

        The observation window is discretised into ``FF_DUR``-sized cells.  Each
        FF scan is assigned to a cell (and a source) such that:

        - exactly ``n_ff`` scans are placed (from ``_determine_ff_count``);
        - each FF lies fully inside one science block (so it cleanly splits it);
        - the FF source is visible by the required number of antennas there;
        - the **minimum gap between consecutive FF scans is maximised** (even
          spread along the observation, ``CP_W_SPREAD``);
        - the **angular separation between each FF and the science target it
          interrupts is minimised** (``CP_W_SEPARATION``);
        - high FF elevation is mildly preferred (``CP_W_ELEVATION``).

        Parameters
        ----------
        ff_names : list[str]
            Candidate FF block names (registered in obs.scans).
        layout : list[ScheduledScanBlock]
            The science-only layout (provides the interrupted-target coords).
        eff_t0 : Time
            Effective observation window start.
        eff_t1 : Time
            Effective observation window end.

        Returns
        -------
        list[tuple[Time, str]]
            (start_time, ff_block_name) pairs in chronological order.  Empty if
            ortools is unavailable or no feasible placement exists.
        """
        if not _HAS_ORTOOLS or not ff_names or not layout:
            return []
        dur = self.FF_DUR
        cell_min = dur.to(u.min).value
        win_min = (eff_t1 - eff_t0).to(u.min).value
        K = int(win_min // cell_min)
        if K < 1:
            return []

        cell_t = [eff_t0 + (c * cell_min) * u.min for c in range(K)]
        min_req = self._n_ant if self.require_all else self.min_ant
        ff_coord = {fn: self._block_source_coord(fn) for fn in ff_names}

        # Candidate (cell, ff-source) placements with separation + elevation cost.
        cand: dict[tuple[int, int], tuple[float, float]] = {}
        cand_by_cell: dict[int, list[int]] = {}
        for c in range(K):
            occ_s = self._occupant_at(layout, cell_t[c])
            occ_e = self._occupant_at(layout, cell_t[c] + dur - 1 * u.s)
            if occ_s is None or occ_e is None or occ_s[0] != occ_e[0]:
                continue  # FF would cross a science-block boundary; skip cell
            occ_coord = occ_s[1]
            for fi, fn in enumerate(ff_names):
                if self._vis_at(fn, cell_t[c]) < min_req:
                    continue
                sep = 0.0
                if occ_coord is not None and ff_coord[fn] is not None:
                    sep = float(ff_coord[fn].separation(occ_coord).deg)
                cand[(c, fi)] = (sep, self._elev_at(fn, cell_t[c]))
                cand_by_cell.setdefault(c, []).append(fi)
        if not cand:
            return []

        usable_cells = sorted(cand_by_cell.keys())
        n_ff = min(self._determine_ff_count((eff_t1 - eff_t0).to(u.h).value), len(usable_cells))
        if n_ff < 1:
            return []

        model = cp_model.CpModel()
        x: dict[tuple[int, int, int], object] = {}
        for i in range(n_ff):
            for (c, fi) in cand:
                x[i, c, fi] = model.NewBoolVar(f"x_{i}_{c}_{fi}")
            model.Add(sum(x[i, c, fi] for (c, fi) in cand) == 1)
        for c in usable_cells:
            model.Add(sum(x[i, c, fi] for i in range(n_ff) for fi in cand_by_cell[c]) <= 1)

        pos = [model.NewIntVar(0, K - 1, f"pos_{i}") for i in range(n_ff)]
        for i in range(n_ff):
            model.Add(pos[i] == sum(c * x[i, c, fi] for (c, fi) in cand))
        for i in range(n_ff - 1):
            model.Add(pos[i + 1] >= pos[i] + 1)

        obj = []
        if n_ff > 1:
            gmin = model.NewIntVar(0, K, "gmin")
            for i in range(n_ff - 1):
                model.Add(gmin <= pos[i + 1] - pos[i])
            obj.append(self.CP_W_SPREAD * gmin)
        else:
            dev = model.NewIntVar(0, K, "dev")
            model.AddAbsEquality(dev, pos[0] - K // 2)
            obj.append(-self.CP_W_SPREAD * dev)

        sep_term = sum(int(round(cand[(c, fi)][0])) * x[i, c, fi]
                       for i in range(n_ff) for (c, fi) in cand)
        elev_term = sum(int(round(cand[(c, fi)][1])) * x[i, c, fi]
                        for i in range(n_ff) for (c, fi) in cand)
        model.Maximize(sum(obj) - self.CP_W_SEPARATION * sep_term + self.CP_W_ELEVATION * elev_term)

        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = self.CP_MAX_SOLVE_SECONDS
        solver.parameters.random_seed = 0
        solver.parameters.num_search_workers = 1
        status = solver.Solve(model)
        if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            log.warning("CP-SAT FF placement found no solution (status %s).", solver.StatusName(status))
            return []

        placements: list[tuple[Time, str]] = []
        for i in range(n_ff):
            for (c, fi) in cand:
                if solver.Value(x[i, c, fi]) == 1:
                    placements.append((cell_t[c], ff_names[fi]))
                    break
        placements.sort(key=lambda p: p[0].mjd)
        log.info("CP-SAT placed %d FF scans (spread=%s, status=%s).",
                 len(placements), 'maximin' if n_ff > 1 else 'centred', solver.StatusName(status))
        return placements

    def _polcal_target_times(self, eff_t0: Time, eff_t1: Time) -> list[tuple[Time, str]]:
        """Create polcal blocks and return their target (time, name) at 10/50/90% of the window.

        Parameters
        ----------
        eff_t0 : Time
        Effective observation window start.
        eff_t1 : Time
        Effective observation window end.

        Returns
        -------
        list[tuple[Time, str]]
        List of (target_time, block_name) tuples.
        """
        names = self._create_polcal_blocks()
        total = (eff_t1 - eff_t0)
        return [(eff_t0 + total * frac, pc) for frac, pc in zip([0.1, 0.5, 0.9], names)]

    def _make_aux_block(self, name: str, a: Time, b: Time, label: str) -> ScheduledScanBlock:
        """Build a single-scan calibrator ScheduledScanBlock (FF/polcal/eMERLIN).

        Parameters
        ----------
        name : str
        Block name in obs.scans.
        a : Time
        Start time.
        b : Time
        End time.
        label : str
        Label for the scheduled block.

        Returns
        -------
        ScheduledScanBlock
        The scheduled calibrator block.
        """
        block = self.obs.scans[name]
        return ScheduledScanBlock(
            name=label, block=block, start_time=a, end_time=b, scans=list(block.scans),
            n_antennas=self._vis_at(name, a), mean_elevation=self._elev_at(name, a))

    def _make_sci_chunk(self, name: str, a: Time, b: Time) -> Optional[ScheduledScanBlock]:
        """Build a science ScheduledScanBlock for [a, b], or None if too short to fit one pass.

        Parameters
        ----------
        name : str
        Block name.
        a : Time
        Start time.
        b : Time
        End time.

        Returns
        -------
        ScheduledScanBlock or None
        Scheduled science block, or None if duration insufficient.
        """
        span = (b - a)
        if span.to(u.min).value <= 0 or span < self._min_block_duration(name):
            return None
        block = self.obs.scans[name]
        try:
            scans = block.fill(span.to(u.min))
        except Exception:
            scans = list(block.scans)
        mid = a + span * 0.5
        return ScheduledScanBlock(
            name=name, block=block, start_time=a, end_time=b, scans=scans,
            n_antennas=self._vis_at(name, mid), mean_elevation=self._elev_at(name, mid))

    def _insert_aux(self, layout: list[ScheduledScanBlock],
                    placements: list[tuple[Time, str, u.Quantity, str]]) -> list[ScheduledScanBlock]:
        """Split science blocks to insert calibrator scans at the requested times.

        Each placement ``(start, block_name, duration, label)`` is assigned to the
        science block that contains its start time and inserted there, splitting
        the surrounding science into the chunks before and after it.  Placements
        that would cross a block boundary or not fit are skipped.

        Parameters
        ----------
        layout : list[ScheduledScanBlock]
            Contiguous science layout.
        placements : list[tuple[Time, str, Quantity, str]]
            Calibrator insertions as (start_time, block_name, duration, label).

        Returns
        -------
        list[ScheduledScanBlock]
            Combined science + calibrator schedule in chronological order.
        """
        if not placements:
            return list(layout)
        buckets: list[list[tuple[Time, str, u.Quantity, str]]] = [[] for _ in layout]
        for p in placements:
            ts = p[0]
            for bi, sb in enumerate(layout):
                if sb.start_time <= ts < sb.end_time:
                    buckets[bi].append(p)
                    break
        counters: dict[str, int] = {}
        result: list[ScheduledScanBlock] = []
        for bi, sb in enumerate(layout):
            plist = sorted(buckets[bi], key=lambda p: p[0].mjd)
            if not plist:
                result.append(sb)
                continue
            cur = sb.start_time
            for (ts, name, pdur, label) in plist:
                ts = max(ts, cur)
                te = ts + pdur
                if te > sb.end_time:
                    continue  # does not fit in the remaining block span
                if ts > cur:
                    chunk = self._make_sci_chunk(sb.name, cur, ts)
                    if chunk:
                        result.append(chunk)
                counters[label] = counters.get(label, 0) + 1
                result.append(self._make_aux_block(name, ts, te, f"{label}_{counters[label]}"))
                cur = te
            if cur < sb.end_time:
                chunk = self._make_sci_chunk(sb.name, cur, sb.end_time)
                if chunk:
                    result.append(chunk)
        return result

    def _schedule_cpsat(self, sci_names: list[str], ff_names: list[str],
                        t0: Time, t1: Time) -> list[ScheduledScanBlock]:
        """Constraint-programming schedule: science layout + CP-SAT FF/calibrator insertion.

        Builds a contiguous science layout, then inserts fringe finders (placed by
        CP-SAT), eMERLIN 3C286, and polcal scans by splitting the science blocks.

        Parameters
        ----------
        sci_names : list[str]
            Science block names.
        ff_names : list[str]
            FF block names.
        t0 : Time
            Observation start time.
        t1 : Time
            Observation end time.

        Returns
        -------
        list[ScheduledScanBlock]
            Full schedule in chronological order, or empty on failure.
        """
        eff_t0, eff_t1 = self._effective_obs_window(t0, t1)
        layout = self._layout_science(sci_names, eff_t0, eff_t1)
        if not layout:
            return []

        aux: list[tuple[Time, str, u.Quantity, str]] = []
        for ts, fn in self._cpsat_ff_placement(ff_names, layout, eff_t0, eff_t1):
            aux.append((ts, fn, self.FF_DUR, 'FF'))
        if self._has_emerlin():
            em = self._create_emerlin_3c286_block()
            if em:
                r = self._find_best(em, eff_t0, eff_t1, self.EMERLIN_3C286_DUR)
                if r:
                    aux.append((r[0], em, self.EMERLIN_3C286_DUR, 'eMERLIN_3C286'))
        if self._polcal:
            for ts, pc in self._polcal_target_times(eff_t0, eff_t1):
                aux.append((ts, pc, self.POLCAL_DUR, 'POLCAL'))

        return sorted(self._insert_aux(layout, aux), key=lambda b: b.start_time.mjd)

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def schedule(self) -> dict[str, ScanBlock]:
        """Generate the complete observation schedule.

        For a single science block: FF (evenly spread) → eMERLIN 3C286 → polcal → science fills gaps.
        For multiple science blocks: multi-source mode with equalised time per source,
        optimal window selection, slew minimisation, and FF scans distributed between
        source blocks.

        Returns
        -------
        dict[str, ScanBlock]
            Ordered mapping ``'NNN_label' → ScanBlock``.
        """
        t0, t1 = self.obs.times[0], self.obs.times[-1]

        sci_names = [n for n in self._sci_blocks()
                     if not n.startswith('FF_') and not n.startswith('POLCAL_') and n != 'eMERLIN_3C286']
        ff_names = self._select_fringefinders()

        # ---- Preferred path: CP-SAT constraint-programming layout ----
        if _HAS_ORTOOLS and sci_names:
            cp_sched = self._schedule_cpsat(sci_names, ff_names, t0, t1)
            if cp_sched:
                self._scheduled = sorted(cp_sched, key=lambda b: b.start_time.mjd)
                return {f"{i + 1:03d}_{b.name}": b.block for i, b in enumerate(self._scheduled)}
            log.warning("CP-SAT scheduling produced no blocks; falling back to heuristic scheduler.")

        if len(sci_names) > 1:
            # ---- Multi-source path ----
            all_sched: list[ScheduledScanBlock] = self._schedule_multi_source(sci_names, ff_names, t0, t1)

            # eMERLIN 3C286 and polcal are added in any remaining gaps
            if self._has_emerlin():
                em_block = self._create_emerlin_3c286_block()
                if em_block:
                    slots = self._get_slots(all_sched, t0, t1)
                    em = self._schedule_single_block(em_block, slots, 'eMERLIN_3C286')
                    if em:
                        all_sched.append(em)
            if self._polcal:
                polcal_names = self._create_polcal_blocks()
                slots = self._get_slots(all_sched, t0, t1)
                all_sched.extend(self._schedule_polcal(polcal_names, slots))
        else:
            # ---- Single-source path (original logic) ----
            ff_sched = self._schedule_ff(ff_names, t0, t1) if ff_names else []
            all_sched = list(ff_sched)

            if self._has_emerlin():
                em_block = self._create_emerlin_3c286_block()
                if em_block:
                    slots = self._get_slots(all_sched, t0, t1)
                    em = self._schedule_single_block(em_block, slots, 'eMERLIN_3C286')
                    if em:
                        all_sched.append(em)

            if self._polcal:
                polcal_names = self._create_polcal_blocks()
                slots = self._get_slots(all_sched, t0, t1)
                all_sched.extend(self._schedule_polcal(polcal_names, slots))

            slots = self._get_slots(all_sched, t0, t1)
            all_sched.extend(self._schedule_sci(sci_names, slots))

        self._scheduled = sorted(all_sched, key=lambda b: b.start_time.mjd)
        return {f"{i + 1:03d}_{b.name}": b.block for i, b in enumerate(self._scheduled)}

    def get_scheduled_blocks(self) -> list[ScheduledScanBlock]:
        """Return scheduled blocks in chronological order.

        Returns
        -------
        list[ScheduledScanBlock]
            Scheduled blocks sorted by start time.
        """
        return self._scheduled

    # ------------------------------------------------------------------
    # Key-file generation  (REWRITTEN — group/rep + Jb1)
    # ------------------------------------------------------------------

    @staticmethod
    def _scans_match(a: list[Scan], b: list[Scan]) -> bool:
        """Return True if two scan lists have the same sources and durations.

        Parameters
        ----------
        a : list[Scan]
        First scan list.
        b : list[Scan]
        Second scan list.

        Returns
        -------
        bool
        True if sources and durations match within 1 second tolerance.
        """
        if len(a) != len(b):
            return False
        return all(x.source.name == y.source.name
                   and abs(x.duration.to(u.s).value - y.duration.to(u.s).value) < 1.0
                   for x, y in zip(a, b))

    @staticmethod
    def _detect_cycle(scans: list[Scan]) -> tuple[list[Scan], int, list[Scan]]:
        """Find the repeating cycle in scans that covers the most total scans.

        Tries all candidate cycle lengths and picks the one where
        ``reps * cycle_len`` is maximised (i.e. covers the largest portion
        of the scan list).  This avoids the greedy trap of matching a short
        sub-cycle (e.g. pcal+target) when a longer cycle (including the
        check source every N) exists.

        Parameters
        ----------
        scans : list[Scan]
        List of scans to analyze.

        Returns
        -------
        tuple[list[Scan], int, list[Scan]]
        (cycle, n_reps, remainder) where cycle is one repetition of the pattern,
        n_reps is the number of full repetitions (≥ 1), and remainder is the
        left-over scans after the last full repetition.
        """
        n = len(scans)
        best_clen, best_reps, best_covered = 0, 1, 0
        for clen in range(2, n // 2 + 1):
            pattern = scans[:clen]
            reps, i = 0, 0
            while i + clen <= n:
                if ObservationScheduler._scans_match(pattern, scans[i:i + clen]):
                    reps += 1
                    i += clen
                else:
                    break
            covered = reps * clen
            if reps >= 2 and covered > best_covered:
                best_clen, best_reps, best_covered = clen, reps, covered
        if best_clen > 0:
            return scans[:best_clen], best_reps, scans[best_covered:]
        return scans, 1, []

    def _jb1_station_name(self) -> Optional[str]:
        """Return the SCHED-style station name for Jb1, or None.

        Returns
        -------
        str or None
        SCHED station name for JB1, or None if not in array.
        """
        for s in self.obs.stations:
            if s.codename.upper() == 'JB1':
                return s.name
        return None

    def _format_scan_line(self, scan: Scan, gap: str = '0:00', indent: str = '',
                          dur_override: Optional[u.Quantity] = None) -> str:
        """Format a single scan as a SCHED source= line.

        Parameters
        ----------
        scan : Scan
        Scan to format.
        gap : str, optional
        Gap before scan in 'M:SS' format. Default is '0:00'.
        indent : str, optional
        Indentation string. Default is ''.
        dur_override : Quantity or None, optional
        If given, use this duration instead of scan.duration.

        Returns
        -------
        str
        Formatted SCHED source line.
        """
        intent = _intent_str(scan.source.type)
        intent_part = f" intent='{intent}'" if intent else ''
        dur = dur_override if dur_override is not None else scan.duration
        return f"{indent}source='{scan.source.name}' gap={gap} dur={_fmt_dur(dur)}{intent_part} /"

    def _stations_line(self, exclude: Optional[set[str]] = None, indent: str = '') -> str:
        """Build a ``stations = ...`` line, optionally excluding codenames.

        Parameters
        ----------
        exclude : set[str] or None, optional
        Set of codenames to exclude. Default is None.
        indent : str, optional
        Indentation string. Default is ''.

        Returns
        -------
        str
        Formatted stations line.
        """
        names = [s.name for s in self.obs.stations
                 if exclude is None or s.codename.upper() not in exclude]
        return f"{indent}stations = {', '.join(names)}"

    def _build_science_scans(self, sb: ScheduledScanBlock) -> list[str]:
        """Convert a science ScheduledScanBlock into SCHED key-file lines.

        Detects repeating cycles and emits ``group N rep R``.
        If Jb1 is in the array, every other PHASECAL scan inside a group
        excludes Jb1 to stay under 12 source-changes / hour.

        Parameters
        ----------
        sb : ScheduledScanBlock
        Scheduled science block to convert.

        Returns
        -------
        list[str]
        List of SCHED key-file lines.
        """
        cycle, n_reps, remainder = self._detect_cycle(sb.scans)
        has_jb1 = self._has_jb1()
        all_stations = self._stations_line()
        lines: list[str] = []

        if n_reps >= 2:
            # Build one cycle with Jb1 handling
            cycle_lines = self._format_cycle(cycle, has_jb1, indent='    ')
            n_scan_lines = sum(1 for ln in cycle_lines if ln.strip().startswith("source="))
            lines.append(f"\ngroup {n_scan_lines} rep {n_reps}")
            lines.append(f"    {all_stations.strip()}")
            lines.extend(cycle_lines)
            # Remainder (partial last cycle) — also needs gaps
            if remainder:
                lines.append(all_stations)
                rem_gaps = self._compute_gap_positions(remainder)
                for idx, scan in enumerate(remainder):
                    gap_str = '0:00'
                    dur_ov: Optional[u.Quantity] = None
                    if idx in rem_gaps:
                        gap_str = _fmt_dur(self.GAP_DUR)
                        dur_ov = scan.duration - self.GAP_DUR
                    lines.append(self._format_scan_line(scan, gap=gap_str, dur_override=dur_ov))
        else:
            # No repeating cycle detected — still insert gaps
            all_gaps = self._compute_gap_positions(sb.scans)
            for idx, scan in enumerate(sb.scans):
                gap_str = '0:00'
                dur_ov = None
                if idx in all_gaps:
                    gap_str = _fmt_dur(self.GAP_DUR)
                    dur_ov = scan.duration - self.GAP_DUR
                lines.append(self._format_scan_line(scan, gap=gap_str, dur_override=dur_ov))
        return lines

    def _compute_gap_positions(self, cycle: list[Scan]) -> set[int]:
        """Determine which scan indices in the cycle should have a gap inserted.

        Distributes gaps evenly across the cycle at approximately GAP_INTERVAL
        spacing.  Only PHASECAL scans whose duration exceeds GAP_DUR are
        eligible.  The gap duration is later subtracted from the scan duration
        so the total time per scan slot stays the same.

        Parameters
        ----------
        cycle : list[Scan]
        List of scans in the cycle.

        Returns
        -------
        set[int]
        Indices into cycle that should get ``gap=0:30``.
        """
        total_dur = sum(s.duration.to(u.min).value for s in cycle)
        interval_min = self.GAP_INTERVAL.to(u.min).value
        gap_dur_s = self.GAP_DUR.to(u.s).value

        if total_dur < interval_min * 0.8:
            return set()

        # Find eligible phasecal positions with their cumulative start times
        pcal_positions: list[tuple[int, float]] = []
        cum = 0.0
        for i, scan in enumerate(cycle):
            if (scan.source.type == SourceType.PHASECAL
                    and scan.duration.to(u.s).value > gap_dur_s):
                pcal_positions.append((i, cum))
            cum += scan.duration.to(u.min).value

        if not pcal_positions:
            return set()

        # Compute how many gaps we need and their ideal positions
        n_gaps = max(1, round(total_dur / interval_min))
        n_gaps = min(n_gaps, len(pcal_positions))
        ideal_times = [(k + 0.5) * total_dur / n_gaps for k in range(n_gaps)]

        # Assign each ideal gap time to the nearest eligible phasecal
        used: set[int] = set()
        gap_positions: set[int] = set()
        for ideal_t in ideal_times:
            best_idx, best_dist = -1, float('inf')
            for scan_idx, cum_t in pcal_positions:
                if scan_idx in used:
                    continue
                dist = abs(cum_t - ideal_t)
                if dist < best_dist:
                    best_dist, best_idx = dist, scan_idx
            if best_idx >= 0:
                gap_positions.add(best_idx)
                used.add(best_idx)
        return gap_positions

    def _format_cycle(self, cycle: list[Scan], handle_jb1: bool, indent: str = '') -> list[str]:
        """Format one cycle of scans for the SCHED key file.

        - Inserts ~30 s gaps on phasecal scans every ~10 min of elapsed cycle
          time.  The gap is subtracted from the scan duration so the total time
          per scan slot stays the same.
        - Every other PHASECAL scan gets a station line excluding Jb1 when
          *handle_jb1* is True (to stay under 12 source-changes / hour).

        Parameters
        ----------
        cycle : list[Scan]
        List of scans in the cycle.
        handle_jb1 : bool
        Whether to handle Jb1 source-change mitigation.
        indent : str, optional
        Indentation string. Default is ''.

        Returns
        -------
        list[str]
        Formatted SCHED key-file lines.
        """
        lines: list[str] = []
        pcal_idx = 0
        all_stations_line = self._stations_line(indent=indent)
        no_jb1_line = self._stations_line(exclude={'JB1'}, indent=indent)
        gap_set = self._compute_gap_positions(cycle)

        for i, scan in enumerate(cycle):
            gap_str = '0:00'
            dur_override: Optional[u.Quantity] = None
            if i in gap_set:
                gap_str = _fmt_dur(self.GAP_DUR)
                dur_override = scan.duration - self.GAP_DUR

            if handle_jb1 and scan.source.type == SourceType.PHASECAL:
                pcal_idx += 1
                if pcal_idx % 2 == 0:
                    lines.append(no_jb1_line)
                    lines.append(self._format_scan_line(scan, gap=gap_str,
                                                        dur_override=dur_override, indent=indent))
                    lines.append(all_stations_line)
                    continue
            lines.append(self._format_scan_line(scan, gap=gap_str,
                                                dur_override=dur_override, indent=indent))
        return lines

    def generate_key_file(self, experiment_code: str = 'EXCODE', pi_name: str = 'PI Name',
                          pi_email: str = 'pi@example.com', pi_institute: str = 'Institute',
                          setup_file: Optional[str] = None, comments: str = '') -> str:
        """Generate a SCHED .key file from the current schedule.

        Uses ``group N rep R`` for repeated science cycles and excludes Jb1
        from every other phase-cal scan when Jb1 is present.

        Parameters
        ----------
        experiment_code : str
        Experiment code for the observation.
        pi_name : str
        Principal investigator name.
        pi_email : str
        Principal investigator email.
        pi_institute : str
        Principal investigator institute.
        setup_file : str or None, optional
        Frequency setup file name. Default is None.
        comments : str, optional
        Additional comments for the key file. Default is ''.

        Returns
        -------
        str
        Complete .key file content.
        """
        with open(resources.files('vlbiplanobs.data').joinpath('key_file.key.template'), 'r') as f:  # type: ignore
            template = f.read()

        # ---- Source catalog ----
        all_sources: dict[str, Source] = {}
        for sb in self._scheduled:
            for scan in sb.scans:
                if scan.source.name not in all_sources:
                    all_sources[scan.source.name] = scan.source
        src_lines = []
        for src in all_sources.values():
            ra = src.coord.ra.to_string(unit=u.hourangle, sep=':', precision=4, pad=True)
            dec = src.coord.dec.to_string(unit=u.degree, sep=':', precision=3, pad=True, alwayssign=True)
            src_lines.append(f"  source='{src.name}' ra={ra} dec={dec} equinox='J2000' /")

        # ---- Scan section ----
        all_stations = self._stations_line()
        scan_lines: list[str] = [all_stations, '']
        for sb in self._scheduled:
            is_ff = sb.block.has(SourceType.FRINGEFINDER)
            is_emerlin = sb.block.has(SourceType.AMPLITUDECAL)
            is_polcal = sb.block.has(SourceType.POLCAL)

            if is_ff or is_emerlin or is_polcal:
                gap = '1:30' if is_emerlin else '0:00'
                for scan in sb.scans:
                    scan_lines.append(self._format_scan_line(scan, gap=gap))
            else:
                scan_lines.extend(self._build_science_scans(sb))
            scan_lines.append('')  # blank line between blocks

        # ---- Template substitution ----
        setup_str = f"setup = '{setup_file}'" if setup_file else "nosetup   ! TODO: Add frequency setup"
        obs_mode = (f"{self.obs.band} {int(self.obs.datarate.to(u.Mbit / u.s).value)} Mbps"
                    if self.obs.band and self.obs.datarate is not None else "VLBI")
        replacements = {
            '{GENERATION_DATE}': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            '{EXPERIMENT_CODE}': experiment_code.upper(),
            '{PI_NAME}': pi_name, '{PI_EMAIL}': pi_email, '{PI_INSTITUTE}': pi_institute,
            '{OBS_MODE}': obs_mode, '{COMMENTS}': comments,
            '{CORAVG}': str(int(self.obs.inttime.to(u.s).value)),
            '{CORCHAN}': str(self.obs.channels) if self.obs.channels else '32',
            '{CORNANT}': str(len(self.obs.stations)),
            '{STATIONS_CATALOG}': 'none',
            '{SOURCES}': '\n'.join(src_lines),
            '{SETUP}': setup_str,
            '{YEAR}': str(self.obs.times[0].datetime.year),
            '{MONTH}': str(self.obs.times[0].datetime.month),
            '{DAY}': str(self.obs.times[0].datetime.day),
            '{START_TIME}': self.obs.times[0].datetime.strftime('%H:%M:%S'),
            '{STATIONS}': ', '.join(s.codename.upper() for s in self.obs.stations),
            '{SCANS}': '\n'.join(scan_lines),
        }
        result = template
        for placeholder, value in replacements.items():
            result = result.replace(placeholder, value)
        return result
