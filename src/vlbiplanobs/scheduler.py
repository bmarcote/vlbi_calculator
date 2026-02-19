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
import logging
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from .sources import Source, SourceType, ScanBlock, Scan
from .observation import Observation

log = logging.getLogger(__name__)

_EMERLIN_CODES = {'CM', 'KN', 'PI', 'DA', 'DE', 'JB2', 'JB1'}
_POLCAL_NAMES = ['3C84', 'OQ208', 'DA193']
_3C286_COORD = '13h31m08.288s +30d30m32.96s'
_JB1_MAX_SRC_CHANGES_PER_HOUR = 12


def _fmt_dur(q: u.Quantity) -> str:
    """Format an astropy duration quantity as 'M:SS' for SCHED key files."""
    total_sec = int(round(q.to(u.s).value))
    return f"{total_sec // 60}:{total_sec % 60:02d}"


def _intent_str(stype: SourceType) -> str:
    """Map a SourceType to a SCHED intent string (empty if none applies)."""
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

    def __init__(self, observation: Observation, min_antennas: int = 2,
                 require_all_antennas: bool = False,
                 fringefinder_spec: Optional[list[str]] = None, polcal: bool = False):
        self.obs = observation
        self.min_ant = min_antennas
        self.require_all = require_all_antennas
        self._ff_spec = fringefinder_spec or ['2']
        self._polcal = polcal
        self._scheduled: list[ScheduledScanBlock] = []
        self._precompute()

    # ------------------------------------------------------------------
    # Pre-computation
    # ------------------------------------------------------------------

    def _precompute(self):
        """Build per-block visibility and mean-elevation arrays over obs.times."""
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
        return [n for n in self._blocks if self._is_ff.get(n, False)]

    def _sci_blocks(self) -> list[str]:
        return [n for n in self._blocks if not self._is_ff.get(n, False)]

    def _t2i(self, t: Time) -> int:
        return int(np.argmin(np.abs(self.obs.times.mjd - t.mjd)))

    def _i2t(self, i: int) -> Time:
        return self.obs.times[min(i, self._n_times - 1)]

    def _vis_at(self, name: str, t: Time) -> int:
        return int(self._vis[name][self._t2i(t)])

    def _elev_at(self, name: str, t: Time) -> float:
        return float(self._elev[name][self._t2i(t)])

    def _vis_mean(self, name: str, t0: Time, t1: Time) -> float:
        """Mean number of antennas visible for *name* between t0 and t1."""
        i0, i1 = self._t2i(t0), self._t2i(t1)
        arr = self._vis[name][i0:max(i1, i0 + 1)]
        return float(np.mean(arr))

    def _elev_mean(self, name: str, t0: Time, t1: Time) -> float:
        """Mean elevation for *name* between t0 and t1."""
        i0, i1 = self._t2i(t0), self._t2i(t1)
        arr = self._elev[name][i0:max(i1, i0 + 1)]
        return float(np.nanmean(arr))

    # ------------------------------------------------------------------
    # Slot helpers
    # ------------------------------------------------------------------

    def _find_best(self, name: str, t0: Time, t1: Time,
                   dur: u.Quantity) -> Optional[tuple[Time, Time, int, float]]:
        """Find optimal placement for *name* in [t0, t1] with given duration.

        Returns (start, end, min_antennas, mean_elevation) or None.
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
        """Return free time slots between already-scheduled blocks."""
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
        """Register a dynamically-added block in pre-computed arrays."""
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

        from . import calibrators

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
        """From a list of CalibratorSource candidates, pick the one closest to science targets."""
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
        """Return the mean sky position of all science targets, or None."""
        ras, decs = [], []
        for name in self._sci_blocks():
            for src in self.obs.scans[name].sources(SourceType.TARGET):
                ras.append(src.coord.ra.deg)
                decs.append(src.coord.dec.deg)
        if not ras:
            return None
        return SkyCoord(float(np.mean(ras)), float(np.mean(decs)), unit=u.deg)

    def _register_ff_sources(self, sources: list[Source]) -> list[str]:
        """Register unique FF Source objects as scan blocks. Returns block names."""
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
        """Look up named fringe-finder sources from the RFC catalog."""
        from . import calibrators
        cat = calibrators.RFCCatalog(min_flux=0.0 * u.Jy, band='c')
        result: list[Source] = []
        for name in names:
            rfc_src = cat.get_source(name)
            if rfc_src is not None:
                result.append(Source(name=rfc_src.name, coordinates=rfc_src.coord,
                                    source_type=SourceType.FRINGEFINDER, other_names=[rfc_src.ivsname]))
            else:
                try:
                    result.append(Source.source_from_str(name, source_type=SourceType.FRINGEFINDER))
                except ValueError:
                    log.warning("Fringe-finder source '%s' not found — skipping.", name)
        return result

    # ------------------------------------------------------------------
    # Polcal / eMERLIN helpers
    # ------------------------------------------------------------------

    def _create_polcal_blocks(self) -> list[str]:
        """Create 3C84 / OQ208 / DA193 polcal blocks. Returns block names."""
        from . import calibrators
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
        return any(s.codename.upper() in _EMERLIN_CODES for s in self.obs.stations)

    def _has_jb1(self) -> bool:
        """Return True if Jodrell Bank Mk2 (Jb1 / JB1) is in the array."""
        return any(s.codename.upper() == 'JB1' for s in self.obs.stations)

    def _create_emerlin_3c286_block(self) -> Optional[str]:
        """Create a 3C286 scan block for eMERLIN flux-scale calibration."""
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
        t0, t1 : Time
            Observation boundaries.

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

        if obs_dur_h <= 1.5:
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
        """Spread polcal scans at 10 %, 50 %, 90 % of available time."""
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
        """Place a single auxiliary block in the best available slot."""
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
    # Main entry point
    # ------------------------------------------------------------------

    def schedule(self) -> dict[str, ScanBlock]:
        """Generate the complete observation schedule.

        Order: FF (evenly spread) → eMERLIN 3C286 → polcal → science fills gaps.

        Returns
        -------
        dict[str, ScanBlock]
            Ordered mapping ``'NNN_label' → ScanBlock``.
        """
        t0, t1 = self.obs.times[0], self.obs.times[-1]

        # 1. Fringe finders — evenly spread, never consecutive
        ff_names = self._select_fringefinders()
        ff_sched = self._schedule_ff(ff_names, t0, t1) if ff_names else []
        all_sched: list[ScheduledScanBlock] = list(ff_sched)

        # 2. eMERLIN 3C286 — place in first available slot after first FF
        if self._has_emerlin():
            em_block = self._create_emerlin_3c286_block()
            if em_block:
                slots = self._get_slots(all_sched, t0, t1)
                em = self._schedule_single_block(em_block, slots, 'eMERLIN_3C286')
                if em:
                    all_sched.append(em)

        # 3. Polcal — spread at 10 / 50 / 90 %
        if self._polcal:
            polcal_names = self._create_polcal_blocks()
            slots = self._get_slots(all_sched, t0, t1)
            all_sched.extend(self._schedule_polcal(polcal_names, slots))

        # 4. Science — fill remaining gaps, optimised per slot
        sci_names = [n for n in self._sci_blocks()
                     if not n.startswith('FF_') and not n.startswith('POLCAL_') and n != 'eMERLIN_3C286']
        slots = self._get_slots(all_sched, t0, t1)
        all_sched.extend(self._schedule_sci(sci_names, slots))

        self._scheduled = sorted(all_sched, key=lambda b: b.start_time.mjd)
        return {f"{i + 1:03d}_{b.name}": b.block for i, b in enumerate(self._scheduled)}

    def get_scheduled_blocks(self) -> list[ScheduledScanBlock]:
        """Return scheduled blocks in chronological order."""
        return self._scheduled

    # ------------------------------------------------------------------
    # Key-file generation  (REWRITTEN — group/rep + Jb1)
    # ------------------------------------------------------------------

    @staticmethod
    def _scans_match(a: list[Scan], b: list[Scan]) -> bool:
        """Return True if two scan lists have the same sources and durations."""
        if len(a) != len(b):
            return False
        return all(x.source.name == y.source.name
                   and abs(x.duration.to(u.s).value - y.duration.to(u.s).value) < 1.0
                   for x, y in zip(a, b))

    @staticmethod
    def _detect_cycle(scans: list[Scan]) -> tuple[list[Scan], int, list[Scan]]:
        """Find the repeating cycle in *scans* that covers the most total scans.

        Tries all candidate cycle lengths and picks the one where
        ``reps * cycle_len`` is maximised (i.e. covers the largest portion
        of the scan list).  This avoids the greedy trap of matching a short
        sub-cycle (e.g. pcal+target) when a longer cycle (including the
        check source every N) exists.

        Returns
        -------
        cycle : list[Scan]
            One repetition of the pattern.
        n_reps : int
            How many full repetitions (≥ 1).
        remainder : list[Scan]
            Left-over scans after the last full repetition.
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
        """Return the SCHED-style station name for Jb1, or None."""
        for s in self.obs.stations:
            if s.codename.upper() == 'JB1':
                return s.name
        return None

    def _format_scan_line(self, scan: Scan, gap: str = '0:00', indent: str = '') -> str:
        """Format a single scan as a SCHED source= line."""
        intent = _intent_str(scan.source.type)
        intent_part = f" intent='{intent}'" if intent else ''
        return f"{indent}source='{scan.source.name}' gap={gap} dur={_fmt_dur(scan.duration)}{intent_part} /"

    def _stations_line(self, exclude: Optional[set[str]] = None, indent: str = '') -> str:
        """Build a ``stations = ...`` line, optionally excluding codenames."""
        names = [s.name for s in self.obs.stations
                 if exclude is None or s.codename.upper() not in exclude]
        return f"{indent}stations = {', '.join(names)}"

    def _build_science_scans(self, sb: ScheduledScanBlock) -> list[str]:
        """Convert a science ScheduledScanBlock into SCHED key-file lines.

        Detects repeating cycles and emits ``group N rep R``.
        If Jb1 is in the array, every other PHASECAL scan inside a group
        excludes Jb1 to stay under 12 source-changes / hour.
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
            # Remainder (partial last cycle)
            if remainder:
                lines.append(all_stations)
                for scan in remainder:
                    lines.append(self._format_scan_line(scan))
        else:
            for scan in sb.scans:
                lines.append(self._format_scan_line(scan))
        return lines

    def _format_cycle(self, cycle: list[Scan], handle_jb1: bool, indent: str = '') -> list[str]:
        """Format one cycle of scans, adding Jb1-exclusion station overrides.

        Every other PHASECAL scan gets a station line excluding Jb1, followed
        by a station line restoring the full array.
        """
        lines: list[str] = []
        pcal_idx = 0
        all_stations_line = self._stations_line(indent=indent)
        no_jb1_line = self._stations_line(exclude={'JB1'}, indent=indent)

        for scan in cycle:
            if handle_jb1 and scan.source.type == SourceType.PHASECAL:
                pcal_idx += 1
                if pcal_idx % 2 == 0:
                    lines.append(no_jb1_line)
                    lines.append(self._format_scan_line(scan, indent=indent))
                    lines.append(all_stations_line)
                    continue
            lines.append(self._format_scan_line(scan, indent=indent))
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
        pi_name, pi_email, pi_institute : str
        setup_file : str or None
        comments : str

        Returns
        -------
        str
            Complete .key file content.
        """
        from importlib import resources
        from datetime import datetime

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
