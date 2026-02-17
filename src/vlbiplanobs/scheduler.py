# -*- coding: utf-8 -*-
# Licensed under GPLv3+ - see LICENSE
"""VLBI Observation Scheduler.

This module provides scheduling functionality for VLBI observations, arranging
scan blocks optimally across the observation time range.

Features
--------
- Fringe finder scheduling: auto-select from calibrators or use user-specified sources.
  2 scans at start, 1 at end, every ~2 h if observation > 3 h.
- Polarization calibrator scheduling: 3C84, OQ208, DA193 spread across observation.
- eMERLIN support: adds a 3C286 scan when eMERLIN antennas are present (if visible).
- Science / target block scheduling: fills remaining time optimally.
- Key file generation for SCHED.

Classes
-------
ScheduledScanBlock
    Dataclass representing a scheduled scan block with timing metadata.
ObservationScheduler
    Main scheduler class for arranging scan blocks.
"""

from dataclasses import dataclass, field
from typing import Optional
import logging
import numpy as np
from astropy import units as u
from astropy.time import Time
from .sources import Source, SourceType, ScanBlock, Scan
from .observation import Observation

log = logging.getLogger(__name__)

# eMERLIN station codenames (upper-cased for comparison)
_EMERLIN_CODES = {'CM', 'KN', 'PI', 'DA', 'DE', 'JB2', 'JB1'}

# Well-known polarisation calibrator names (looked up in RFC)
_POLCAL_NAMES = ['3C84', 'OQ208', 'DA193']

# 3C286 coordinates (J2000) – used for eMERLIN flux-scale calibration
_3C286_COORD = '13h31m08.288s +30d30m32.96s'


@dataclass
class ScheduledScanBlock:
    """A scan block scheduled at a specific time with timing and quality metrics.

    Attributes
    ----------
    name : str
        Identifier for this scheduled instance (e.g., 'FF_start_1', 'Target_A').
    block : ScanBlock
        The original ScanBlock being scheduled.
    start_time : Time
        Start time of this scheduled block.
    end_time : Time
        End time of this scheduled block.
    scans : list[Scan]
        The expanded list of scans filling this time slot.
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


class ObservationScheduler:
    """Schedules VLBI scan blocks across the observation time range.

    Arranges fringe finders, polarisation calibrators, eMERLIN flux calibrator,
    and science blocks following standard VLBI conventions, optimising for
    antenna participation and source elevation.

    Parameters
    ----------
    observation : Observation
        The observation containing scan blocks to schedule.
    min_antennas : int, default=2
        Minimum antennas required for valid observation.
    require_all_antennas : bool, default=False
        If True, only schedule when all antennas can observe.
    fringefinder_spec : list[str] or None
        User specification for fringe finders: either source names or a single
        number string (e.g. ['2']) meaning how many FF scans to create.
    polcal : bool, default=False
        Whether polarisation calibration scans are required.
    """

    FF_DUR = 5 * u.min
    FF_INTERVAL = 2 * u.h
    POLCAL_DUR = 5 * u.min
    EMERLIN_3C286_DUR = 5 * u.min

    def __init__(self, observation: Observation, min_antennas: int = 2,
                 require_all_antennas: bool = False,
                 fringefinder_spec: Optional[list[str]] = None,
                 polcal: bool = False):
        self.obs = observation
        self.min_ant = min_antennas
        self.require_all = require_all_antennas
        self._ff_spec = fringefinder_spec or ['2']
        self._polcal = polcal
        self._scheduled: list[ScheduledScanBlock] = []
        self._precompute()

    # ------------------------------------------------------------------
    # Pre-computation helpers
    # ------------------------------------------------------------------

    def _precompute(self):
        """Pre-compute visibility and elevation arrays for all scan blocks.

        Populates internal dictionaries ``_vis``, ``_elev``, and ``_is_ff``
        with time-series data for each block.
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
                               if name in is_obs
                               else np.zeros(self._n_times, dtype=int))
            targets = (block.sources(SourceType.TARGET)
                       or block.sources(SourceType.FRINGEFINDER)
                       or block.sources())
            if targets and targets[0].name in elevs:
                self._elev[name] = np.nanmean(
                    [elevs[targets[0].name][a].value
                     for a in elevs[targets[0].name]], axis=0)
            else:
                self._elev[name] = np.zeros(self._n_times)

    def _ff_blocks(self) -> list[str]:
        """Return names of scan blocks that contain fringe-finder sources."""
        return [n for n in self._blocks if self._is_ff[n]]

    def _sci_blocks(self) -> list[str]:
        """Return names of scan blocks that are science (non-FF) blocks."""
        return [n for n in self._blocks if not self._is_ff[n]]

    def _t2i(self, t: Time) -> int:
        """Convert an astropy Time to the nearest index in obs.times."""
        return int(np.argmin(np.abs(self.obs.times.mjd - t.mjd)))

    def _i2t(self, i: int) -> Time:
        """Convert a time index to an astropy Time."""
        return self.obs.times[min(i, self._n_times - 1)]

    def _vis_at(self, name: str, t: Time) -> int:
        """Number of antennas that can observe *name* at time *t*."""
        return int(self._vis[name][self._t2i(t)])

    def _elev_at(self, name: str, t: Time) -> float:
        """Mean elevation of *name* at time *t*."""
        return float(self._elev[name][self._t2i(t)])

    # ------------------------------------------------------------------
    # Slot / placement helpers
    # ------------------------------------------------------------------

    def _find_best(self, name: str, t0: Time, t1: Time,
                   dur: u.Quantity) -> Optional[tuple[Time, Time, int, float]]:
        """Find optimal time slot for a block within a time range.

        Parameters
        ----------
        name : str
            Block name to schedule.
        t0 : Time
            Earliest allowed start time.
        t1 : Time
            Latest allowed end time.
        dur : Quantity
            Required duration for the block.

        Returns
        -------
        tuple or None
            (start_time, end_time, n_antennas, mean_elevation) if found, else None.
        """
        i0 = self._t2i(t0)
        i1 = self._t2i(t1)
        steps = max(1, int(np.ceil(dur.to(u.min).value / self._dt)))
        if i1 - i0 < steps:
            return None

        vis = self._vis[name]
        elev = self._elev[name]
        min_req = self._n_ant if self.require_all else self.min_ant
        best_score, best_i = -1, -1
        for i in range(i0, i1 - steps + 1):
            n = int(np.min(vis[i:i + steps]))
            e = float(np.mean(elev[i:i + steps]))
            if n >= min_req and not np.isnan(e):
                s = int(n * 100 + e)
                if s > best_score:
                    best_score, best_i = s, i
        if best_i < 0:
            return None

        return (self._i2t(best_i), self._i2t(best_i + steps),
                int(np.min(vis[best_i:best_i + steps])),
                float(np.mean(elev[best_i:best_i + steps])))

    def _get_slots(self, sched: list[ScheduledScanBlock],
                   t0: Time, t1: Time) -> list[tuple[Time, Time]]:
        """Get available time slots between already scheduled blocks.

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
            Available (start, end) time slots.
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

    # ------------------------------------------------------------------
    # Fringe-finder resolution
    # ------------------------------------------------------------------

    def _resolve_ff_blocks(self) -> list[str]:
        """Resolve fringe-finder specification into scan-block names.

        If the observation already has FF blocks in its scans dict, use those.
        Otherwise, auto-select from the calibrators module based on
        ``_ff_spec`` (either source names or a count).

        Returns
        -------
        list[str]
            Names of FF scan blocks present in ``self.obs.scans``.
        """
        existing_ff = self._ff_blocks()
        if existing_ff:
            return existing_ff

        # Determine how many FF scans / which sources to use
        spec = self._ff_spec
        if len(spec) == 1 and spec[0].isdigit():
            n_ff = int(spec[0])
            ff_sources = self._auto_select_fringefinders(n_ff)
        else:
            ff_sources = self._lookup_ff_sources(spec)

        # Register as scan blocks in the observation
        for src in ff_sources:
            block_name = f"FF_{src.name}"
            self.obs.scans[block_name] = ScanBlock(
                [Scan(src, duration=self.FF_DUR)])
            self._register_block(block_name)

        return self._ff_blocks()

    def _auto_select_fringefinders(self, n_scans: int) -> list[Source]:
        """Auto-select fringe-finder sources visible during the observation.

        Uses the calibrators module.  If a single source is visible for the
        entire observation it is reused; otherwise different sources may be
        selected for different time segments.

        Parameters
        ----------
        n_scans : int
            How many fringe-finder scan slots to fill.

        Returns
        -------
        list[Source]
            Unique fringe-finder Source objects.
        """
        from . import calibrators

        sources_found, _, _ = calibrators.get_fringe_finder_sources(
            self.obs.stations, self.obs.times,
            min_elevation=20 * u.deg, min_flux=1.0 * u.Jy,
            require_all_stations=True)

        if not sources_found:
            sources_found, _, _ = calibrators.get_fringe_finder_sources(
                self.obs.stations, self.obs.times,
                min_elevation=20 * u.deg, min_flux=0.5 * u.Jy,
                require_all_stations=False)

        if not sources_found:
            log.warning("No fringe-finder sources found for this observation.")
            return []

        # If the top source is visible all the time, reuse it
        result: list[Source] = []
        seen_names: set[str] = set()
        for src in sources_found:
            if len(result) >= n_scans:
                break
            if src.name not in seen_names:
                ff_src = Source(name=src.name, coordinates=src.coord,
                               source_type=SourceType.FRINGEFINDER,
                               other_names=[src.ivsname])
                result.append(ff_src)
                seen_names.add(src.name)

        # If we still don't have enough, reuse the first one
        while len(result) < n_scans and result:
            result.append(result[0])

        return result

    def _lookup_ff_sources(self, names: list[str]) -> list[Source]:
        """Look up named fringe-finder sources from the RFC catalog.

        Parameters
        ----------
        names : list[str]
            Source names (J2000 or IVS).

        Returns
        -------
        list[Source]
            Resolved Source objects with type FRINGEFINDER.
        """
        from . import calibrators

        cat = calibrators.RFCCatalog(min_flux=0.0 * u.Jy, band='c')
        result: list[Source] = []
        for name in names:
            rfc_src = cat.get_source(name)
            if rfc_src is not None:
                result.append(Source(name=rfc_src.name, coordinates=rfc_src.coord,
                                    source_type=SourceType.FRINGEFINDER,
                                    other_names=[rfc_src.ivsname]))
            else:
                try:
                    result.append(Source.source_from_str(name, source_type=SourceType.FRINGEFINDER))
                except ValueError:
                    log.warning("Fringe-finder source '%s' not found — skipping.", name)
        return result

    def _register_block(self, name: str):
        """Register a newly-added scan block in the pre-computed arrays.

        Parameters
        ----------
        name : str
            Name of the block already added to ``self.obs.scans``.
        """
        if name in self._vis:
            return
        self._blocks.append(name)
        block = self.obs.scans[name]
        self._is_ff[name] = block.has(SourceType.FRINGEFINDER)
        # Approximate visibility: assume visible everywhere (will be refined at placement)
        self._vis[name] = np.full(self._n_times, self._n_ant, dtype=int)
        self._elev[name] = np.full(self._n_times, 45.0)

    # ------------------------------------------------------------------
    # Polcal / eMERLIN helpers
    # ------------------------------------------------------------------

    def _create_polcal_blocks(self) -> list[str]:
        """Create polarisation calibrator scan blocks (3C84, OQ208, DA193).

        Returns
        -------
        list[str]
            Block names added to the observation.
        """
        from . import calibrators

        cat = calibrators.RFCCatalog(min_flux=0.0 * u.Jy, band='c')
        added: list[str] = []
        for polcal_name in _POLCAL_NAMES:
            rfc_src = cat.get_source(polcal_name)
            if rfc_src is None:
                continue
            src = Source(name=rfc_src.name, coordinates=rfc_src.coord,
                         source_type=SourceType.POLCAL,
                         other_names=[rfc_src.ivsname])
            block_name = f"POLCAL_{rfc_src.ivsname}"
            self.obs.scans[block_name] = ScanBlock(
                [Scan(src, duration=self.POLCAL_DUR)])
            self._register_block(block_name)
            added.append(block_name)
        return added

    def _has_emerlin(self) -> bool:
        """Check if any eMERLIN antennas are in the observation."""
        for s in self.obs.stations:
            if s.codename.upper() in _EMERLIN_CODES:
                return True
        return False

    def _create_emerlin_3c286_block(self) -> Optional[str]:
        """Create a 3C286 scan block for eMERLIN flux-scale calibration.

        Returns
        -------
        str or None
            Block name, or None if 3C286 is not visible.
        """
        src_3c286 = Source(name='3C286', coordinates=_3C286_COORD,
                           source_type=SourceType.AMPLITUDECAL)
        block_name = 'eMERLIN_3C286'
        self.obs.scans[block_name] = ScanBlock(
            [Scan(src_3c286, duration=self.EMERLIN_3C286_DUR)])
        self._register_block(block_name)
        return block_name

    # ------------------------------------------------------------------
    # Scheduling algorithms
    # ------------------------------------------------------------------

    def _schedule_ff(self, ff_names: list[str],
                     t0: Time, t1: Time) -> list[ScheduledScanBlock]:
        """Schedule fringe-finder blocks following VLBI conventions.

        Convention: 2 scans at start, 1 at end, additional every ~2 h for
        observations longer than 3 h.

        Parameters
        ----------
        ff_names : list[str]
            Names of fringe-finder scan blocks.
        t0 : Time
            Observation start time.
        t1 : Time
            Observation end time.

        Returns
        -------
        list[ScheduledScanBlock]
            Scheduled FF blocks.
        """
        if not ff_names:
            return []

        dur = self.FF_DUR
        # Build list of (start, end, label) tuples
        time_slots: list[tuple[Time, Time, str]] = []

        # 2 FF scans at the very beginning
        time_slots.append((t0, t0 + dur, "FF_start_1"))
        time_slots.append((t0 + dur, t0 + 2 * dur, "FF_start_2"))

        # Intermediate FF scans every ~2 h for long observations
        if (t1 - t0).to(u.h) > 3 * u.h:
            t_cur = t0 + 2 * dur
            remaining = t1 - t0 - 3 * dur
            n_mid = max(0, int(remaining.to(u.h).value /
                               self.FF_INTERVAL.to(u.h).value))
            if n_mid > 0:
                gap = remaining / (n_mid + 1)
                for i in range(n_mid):
                    ts = t_cur + gap * (i + 1)
                    time_slots.append((ts, ts + dur, f"FF_mid_{i + 1}"))

        # 1 FF scan at the end
        time_slots.append((t1 - dur, t1, "FF_end"))

        # Distribute across available FF blocks (round-robin if multiple)
        result: list[ScheduledScanBlock] = []
        for idx, (ts, te, label) in enumerate(time_slots):
            block_name = ff_names[idx % len(ff_names)]
            block = self.obs.scans[block_name]
            # Use scans directly for calibrator blocks (avoids float precision in fill())
            result.append(ScheduledScanBlock(
                name=label, block=block, start_time=ts, end_time=te,
                scans=list(block.scans),
                n_antennas=self._vis_at(block_name, ts),
                mean_elevation=self._elev_at(block_name, ts)))
        return result

    def _schedule_polcal(self, polcal_names: list[str],
                         slots: list[tuple[Time, Time]]) -> list[ScheduledScanBlock]:
        """Schedule polcal scans spread as far apart as possible.

        Parameters
        ----------
        polcal_names : list[str]
            Names of polcal scan blocks.
        slots : list[tuple[Time, Time]]
            Available time slots.

        Returns
        -------
        list[ScheduledScanBlock]
            Scheduled polcal blocks.
        """
        if not polcal_names or not slots:
            return []

        # Find 2–3 well-separated placements within available slots
        n_polcal = min(3, len(polcal_names))
        dur = self.POLCAL_DUR
        placed: list[ScheduledScanBlock] = []

        # Spread polcal scans across the observation: start, middle, end
        total_available = sum(((s1 - s0).to(u.min).value for s0, s1 in slots)) * u.min
        if total_available < dur:
            return []

        targets = [0.1, 0.5, 0.9][:n_polcal]
        for frac, pc_name in zip(targets, polcal_names):
            block = self.obs.scans[pc_name]
            # Find the slot that overlaps the desired fractional time
            elapsed = 0.0 * u.min
            target_time_min = frac * total_available
            for s0, s1 in slots:
                slot_dur = (s1 - s0).to(u.min)
                if elapsed + slot_dur >= target_time_min and slot_dur >= dur:
                    offset = max(target_time_min - elapsed, 0.0 * u.min)
                    ts = s0 + offset
                    te = ts + dur
                    if te <= s1:
                        placed.append(ScheduledScanBlock(
                            name=f"POLCAL_{pc_name.split('_')[-1]}",
                            block=block, start_time=ts, end_time=te,
                            scans=[block.scans[0]],
                            n_antennas=self._vis_at(pc_name, ts),
                            mean_elevation=self._elev_at(pc_name, ts)))
                        break
                elapsed += slot_dur
        return placed

    def _schedule_single_block(self, block_name: str,
                               slots: list[tuple[Time, Time]],
                               label: str) -> Optional[ScheduledScanBlock]:
        """Schedule a single auxiliary block (e.g. eMERLIN 3C286) in the best slot.

        Parameters
        ----------
        block_name : str
            Name of the scan block.
        slots : list[tuple[Time, Time]]
            Available time slots.
        label : str
            Human-readable label for the scheduled block.

        Returns
        -------
        ScheduledScanBlock or None
            The scheduled block, or None if it couldn't be placed.
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
            return ScheduledScanBlock(
                name=label, block=block, start_time=ts, end_time=te,
                scans=[block.scans[0]], n_antennas=n, mean_elevation=e)
        return None

    def _schedule_sci(self, names: list[str],
                      slots: list[tuple[Time, Time]]) -> list[ScheduledScanBlock]:
        """Schedule science blocks in available time slots.

        Optimises placement for maximum antenna participation and elevation.

        Parameters
        ----------
        names : list[str]
            Names of science blocks to schedule.
        slots : list[tuple[Time, Time]]
            Available time slots.

        Returns
        -------
        list[ScheduledScanBlock]
            Scheduled science blocks.
        """
        scheduled: list[ScheduledScanBlock] = []
        remaining = list(slots)
        for name in names:
            block = self.obs.scans[name]
            min_dur = sum(s.duration.to(u.min).value for s in block.scans) * u.min
            best = None
            for s0, s1 in remaining:
                if (s1 - s0).to(u.min) >= min_dur:
                    r = self._find_best(name, s0, s1, (s1 - s0).to(u.min))
                    if r and (not best or r[2] * 100 + r[3] > best[2] * 100 + best[3]):
                        best = r
            if best:
                t0, t1, n, e = best
                scheduled.append(ScheduledScanBlock(
                    name=name, block=block, start_time=t0, end_time=t1,
                    scans=block.fill((t1 - t0).to(u.min)),
                    n_antennas=n, mean_elevation=e))
                remaining = (
                    [(a, b) for a, b in remaining if b <= t0 or a >= t1] +
                    [(a, t0) for a, b in remaining if a < t0 < b] +
                    [(t1, b) for a, b in remaining if a < t1 < b])
        return scheduled

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def schedule(self) -> dict[str, ScanBlock]:
        """Generate the complete observation schedule.

        Orchestrates fringe-finder, polcal, eMERLIN 3C286, and science block
        scheduling.  Calibrator scans are spread as far apart as possible,
        science blocks fill the remaining time.

        Returns
        -------
        dict[str, ScanBlock]
            Ordered dictionary with keys like '001_FF_start_1' and ScanBlock
            values, covering the entire observation time range.
        """
        t0, t1 = self.obs.times[0], self.obs.times[-1]

        # 1. Resolve & schedule fringe finders
        ff_names = self._resolve_ff_blocks()
        ff_sched = self._schedule_ff(ff_names, t0, t1)
        all_sched: list[ScheduledScanBlock] = list(ff_sched)

        # 2. Polcal blocks (spread across observation)
        if self._polcal:
            polcal_names = self._create_polcal_blocks()
            slots = self._get_slots(all_sched, t0, t1)
            polcal_sched = self._schedule_polcal(polcal_names, slots)
            all_sched.extend(polcal_sched)

        # 3. eMERLIN 3C286
        if self._has_emerlin():
            em_block = self._create_emerlin_3c286_block()
            if em_block:
                slots = self._get_slots(all_sched, t0, t1)
                em_sched = self._schedule_single_block(em_block, slots, 'eMERLIN_3C286')
                if em_sched:
                    all_sched.append(em_sched)

        # 4. Science / target blocks fill remaining time
        sci_names = self._sci_blocks()
        # Filter out blocks already scheduled (FF, polcal, etc.)
        already = {b.name for b in all_sched}
        sci_names = [n for n in sci_names
                     if n not in already and not n.startswith('FF_')
                     and not n.startswith('POLCAL_') and n != 'eMERLIN_3C286']
        slots = self._get_slots(all_sched, t0, t1)
        sci_sched = self._schedule_sci(sci_names, slots)
        all_sched.extend(sci_sched)

        self._scheduled = sorted(all_sched, key=lambda b: b.start_time.mjd)
        return {f"{i + 1:03d}_{b.name}": b.block
                for i, b in enumerate(self._scheduled)}

    def get_scheduled_blocks(self) -> list[ScheduledScanBlock]:
        """Return the list of scheduled blocks with full timing metadata.

        Returns
        -------
        list[ScheduledScanBlock]
            All scheduled blocks in chronological order.
        """
        return self._scheduled

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def print_schedule(self):
        """Print the schedule in a human-readable format to stdout."""
        s = self.get_scheduled_blocks()
        if not s:
            print("No schedule. Run schedule() first.")
            return
        print("\n" + "=" * 70)
        print("VLBI OBSERVATION SCHEDULE")
        print("=" * 70)
        for i, b in enumerate(s):
            src_names = ', '.join(sc.source.name for sc in b.scans) if b.scans else '—'
            print(f"\n{i + 1:3d}. {b.name}")
            print(f"     {b.start_time.iso}  –  {b.end_time.iso}")
            print(f"     {b.duration.value:.1f} min | {b.n_antennas} ant "
                  f"| {b.mean_elevation:.1f}° | {len(b.scans)} scans")
            print(f"     Sources: {src_names}")
        print("\n" + "=" * 70)

    def generate_key_file(self, experiment_code: str = 'EXCODE',
                          pi_name: str = 'PI Name',
                          pi_email: str = 'pi@example.com',
                          pi_institute: str = 'Institute',
                          setup_file: Optional[str] = None,
                          comments: str = '') -> str:
        """Generate a SCHED .key file from the current schedule.

        Parameters
        ----------
        experiment_code : str
            The experiment code (e.g. 'EG123A').
        pi_name, pi_email, pi_institute : str
            PI metadata.
        setup_file : str or None
            Frequency setup filename.
        comments : str
            Extra comments for the cover letter.

        Returns
        -------
        str
            Complete .key file content.
        """
        from importlib import resources
        from datetime import datetime

        template_path = resources.files('vlbiplanobs.data').joinpath(
            'key_file.key.template')
        with open(template_path, 'r') as f:  # type: ignore
            template = f.read()

        # Collect all unique sources from scheduled blocks
        all_sources: dict[str, Source] = {}
        for sb in self._scheduled:
            for scan in sb.scans:
                if scan.source.name not in all_sources:
                    all_sources[scan.source.name] = scan.source

        sources_lines = []
        for src in all_sources.values():
            ra_str = src.coord.ra.to_string(
                unit=u.hourangle, sep=':', precision=4, pad=True)
            dec_str = src.coord.dec.to_string(
                unit=u.degree, sep=':', precision=3, pad=True, alwayssign=True)
            sources_lines.append(
                f"  source='{src.name}' ra={ra_str} dec={dec_str} "
                f"equinox='J2000' /")
        sources_str = '\n'.join(sources_lines)

        stations_list = ', '.join(
            s.codename.upper() for s in self.obs.stations)

        scans_lines = []
        for sb in self._scheduled:
            for scan in sb.scans:
                dur_min = int(scan.duration.to(u.min).value)
                dur_sec = int(scan.duration.to(u.s).value % 60)
                scans_lines.append(
                    f"source='{scan.source.name}' gap=0:00 "
                    f"dur={dur_min}:{dur_sec:02d} /")
        scans_str = '\n'.join(scans_lines)

        setup_str = (f"setup = '{setup_file}'" if setup_file
                     else "nosetup   ! TODO: Add frequency setup")
        obs_mode = (f"{self.obs.band} "
                    f"{int(self.obs.datarate.to(u.Mbit / u.s).value)} Mbps"
                    if self.obs.band and self.obs.datarate is not None
                    else "VLBI")

        replacements = {
            '{GENERATION_DATE}': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            '{EXPERIMENT_CODE}': experiment_code.upper(),
            '{PI_NAME}': pi_name,
            '{PI_EMAIL}': pi_email,
            '{PI_INSTITUTE}': pi_institute,
            '{OBS_MODE}': obs_mode,
            '{COMMENTS}': comments,
            '{CORAVG}': str(int(self.obs.inttime.to(u.s).value)),
            '{CORCHAN}': str(self.obs.channels) if self.obs.channels else '32',
            '{CORNANT}': str(len(self.obs.stations)),
            '{STATIONS_CATALOG}': 'none',
            '{SOURCES}': sources_str,
            '{SETUP}': setup_str,
            '{YEAR}': str(self.obs.times[0].datetime.year),
            '{MONTH}': str(self.obs.times[0].datetime.month),
            '{DAY}': str(self.obs.times[0].datetime.day),
            '{START_TIME}': self.obs.times[0].datetime.strftime('%H:%M:%S'),
            '{STATIONS}': stations_list,
            '{SCANS}': scans_str,
        }

        key_content = template
        for placeholder, value in replacements.items():
            key_content = key_content.replace(placeholder, value)
        return key_content
