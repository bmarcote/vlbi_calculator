# -*- coding: utf-8 -*-
# Licensed under GPLv3+ - see LICENSE
"""VLBI Observation Scheduler.

This module provides scheduling functionality for VLBI observations, arranging
scan blocks optimally across the observation time range.

Features
--------
- Fringe finder scheduling: 2 at start, 1 at end, every ~2h if observation > 3h
- Science block scheduling: optimized for antenna participation and elevation
- Support for all-antenna constraints

Classes
-------
ScheduledScanBlock
    Dataclass representing a scheduled scan block with timing metadata.
ObservationScheduler
    Main scheduler class for arranging scan blocks.
"""

from dataclasses import dataclass, field
from typing import Optional
import numpy as np
from astropy import units as u
from astropy.time import Time
from .sources import SourceType, ScanBlock, Scan
from .observation import Observation


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
        return (self.end_time - self.start_time).to(u.min)


class ObservationScheduler:
    """Schedules VLBI scan blocks across the observation time range.

    Arranges fringe finders and science blocks following standard VLBI conventions,
    optimizing for antenna participation and source elevation.

    Parameters
    ----------
    observation : Observation
        The observation containing scan blocks to schedule.
    min_antennas : int, default=2
        Minimum antennas required for valid observation.
    require_all_antennas : bool, default=False
        If True, only schedule when all antennas can observe.

    Attributes
    ----------
    FF_DUR : Quantity
        Fringe finder scan duration (5 minutes).
    FF_INTERVAL : Quantity
        Interval between fringe finders for long observations (2 hours).

    Examples
    --------
    >>> scheduler = ObservationScheduler(observation, min_antennas=3)
    >>> schedule = scheduler.schedule()
    >>> scheduler.print_schedule()
    """

    FF_DUR, FF_INTERVAL = 5 * u.min, 2 * u.h

    def __init__(self, observation: Observation, min_antennas: int = 2, require_all_antennas: bool = False):
        self.obs, self.min_ant, self.require_all = observation, min_antennas, require_all_antennas
        self._precompute()

    def _precompute(self):
        """Pre-compute visibility and elevation arrays for all scan blocks.

        Populates internal dictionaries `_vis`, `_elev`, and `_is_ff` with
        time-series data for each block.
        """
        self._blocks = list(self.obs.scans.keys())
        self._n_times, self._n_ant = len(self.obs.times), len(self.obs.stations)
        self._dt = (self.obs.times[1] - self.obs.times[0]).to(u.min).value
        is_obs, elevs = self.obs.is_observable(), self.obs.elevations()
        self._vis, self._elev, self._is_ff = {}, {}, {}

        for name in self._blocks:
            block = self.obs.scans[name]
            self._is_ff[name] = block.has(SourceType.FRINGEFINDER)
            self._vis[name] = np.sum(np.array(list(is_obs[name].values())), axis=0) if name in is_obs else np.zeros(self._n_times, dtype=int)
            targets = block.sources(SourceType.TARGET) or block.sources(SourceType.FRINGEFINDER) or block.sources()
            if targets and targets[0].name in elevs:
                self._elev[name] = np.nanmean([elevs[targets[0].name][a].value for a in elevs[targets[0].name]], axis=0)
            else:
                self._elev[name] = np.zeros(self._n_times)

    def _ff_blocks(self) -> list[str]:
        return [n for n in self._blocks if self._is_ff[n]]

    def _sci_blocks(self) -> list[str]:
        return [n for n in self._blocks if not self._is_ff[n]]

    def _t2i(self, t: Time) -> int:
        return int(np.argmin(np.abs(self.obs.times.mjd - t.mjd)))

    def _i2t(self, i: int) -> Time:
        return self.obs.times[min(i, self._n_times - 1)]

    def _vis_at(self, name: str, t: Time) -> int:
        return int(self._vis[name][self._t2i(t)])

    def _elev_at(self, name: str, t: Time) -> float:
        return float(self._elev[name][self._t2i(t)])

    def _find_best(self, name: str, t0: Time, t1: Time, dur: u.Quantity) -> Optional[tuple[Time, Time, int, float]]:
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
        i0, i1, steps = self._t2i(t0), self._t2i(t1), max(1, int(np.ceil(dur.to(u.min).value / self._dt)))
        if i1 - i0 < steps:
            return None

        vis, elev, min_req = self._vis[name], self._elev[name], self._n_ant if self.require_all else self.min_ant
        best_score, best_i = -1, -1
        for i in range(i0, i1 - steps + 1):
            n, e = int(np.min(vis[i:i+steps])), float(np.mean(elev[i:i+steps]))
            if n >= min_req and not np.isnan(e) and (s := int(n * 100 + e)) > best_score:
                best_score, best_i = s, i
        if best_i < 0:
            return None

        return self._i2t(best_i), self._i2t(best_i + steps), int(np.min(vis[best_i:best_i+steps])), float(np.mean(elev[best_i:best_i+steps]))

    def _schedule_ff(self, ff_names: list[str], t0: Time, t1: Time) -> list[ScheduledScanBlock]:
        """Schedule fringe finder blocks according to VLBI conventions.

        Parameters
        ----------
        ff_names : list[str]
            Names of fringe finder blocks.
        t0 : Time
            Observation start time.
        t1 : Time
            Observation end time.

        Returns
        -------
        list[ScheduledScanBlock]
            Scheduled fringe finder blocks (2 at start, 1 at end, ~2h intervals if >3h).
        """
        if not ff_names:
            return []
        block, dur = self.obs.scans[ff_names[0]], self.FF_DUR
        times = [(t0 + i*dur, t0 + (i+1)*dur, f"FF_start_{i+1}") for i in range(2)]
        if (t1 - t0).to(u.h) > 3 * u.h:
            t_cur, remaining = t0 + 2*dur, t1 - t0 - 3*dur
            n_mid = max(0, int(remaining.to(u.h).value / self.FF_INTERVAL.to(u.h).value))
            if n_mid > 0:
                gap = remaining / (n_mid + 1)
                times += [(t_cur + gap*(i+1), t_cur + gap*(i+1) + dur, f"FF_mid_{i+1}") for i in range(n_mid)]
        times.append((t1 - dur, t1, "FF_end"))
        return [ScheduledScanBlock(name=n, block=block, start_time=s, end_time=e,
                scans=block.fill((e-s).to(u.min)), n_antennas=self._vis_at(ff_names[0], s),
                mean_elevation=self._elev_at(ff_names[0], s)) for s, e, n in times]

    def _get_slots(self, sched: list[ScheduledScanBlock], t0: Time, t1: Time) -> list[tuple[Time, Time]]:
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
        slots = ([(t0, s[0].start_time)] if s[0].start_time > t0 else []) + \
                [(s[i].end_time, s[i+1].start_time) for i in range(len(s)-1) if s[i+1].start_time > s[i].end_time] + \
                ([(s[-1].end_time, t1)] if s[-1].end_time < t1 else [])
        return slots

    def _schedule_sci(self, names: list[str], slots: list[tuple[Time, Time]]) -> list[ScheduledScanBlock]:
        """Schedule science blocks in available time slots.

        Optimizes placement for maximum antenna participation and elevation.

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
        scheduled, remaining = [], list(slots)
        for name in names:
            block = self.obs.scans[name]
            min_dur = sum(s.duration.to(u.min).value for s in block.scans) * u.min
            best = None
            for s0, s1 in remaining:
                if (s1 - s0).to(u.min) >= min_dur and (r := self._find_best(name, s0, s1, (s1-s0).to(u.min))):
                    if not best or r[2]*100 + r[3] > best[2]*100 + best[3]:
                        best = r
            if best:
                t0, t1, n, e = best
                scheduled.append(ScheduledScanBlock(name=name, block=block, start_time=t0, end_time=t1,
                                  scans=block.fill((t1-t0).to(u.min)), n_antennas=n, mean_elevation=e))
                remaining = [(a, b) for a, b in remaining if b <= t0 or a >= t1] + \
                            [(a, t0) for a, b in remaining if a < t0 < b] + \
                            [(t1, b) for a, b in remaining if a < t1 < b]
        return scheduled

    def schedule(self) -> dict[str, ScanBlock]:
        """Generate the complete observation schedule.

        Returns
        -------
        dict[str, ScanBlock]
            Ordered dictionary with keys like '001_FF_start_1' and ScanBlock values,
            covering the entire observation time range.
        """
        t0, t1 = self.obs.times[0], self.obs.times[-1]
        ff = self._schedule_ff(self._ff_blocks(), t0, t1)
        sci = self._schedule_sci(self._sci_blocks(), self._get_slots(ff, t0, t1))
        self._scheduled = sorted(ff + sci, key=lambda b: b.start_time.mjd)
        return {f"{i+1:03d}_{b.name}": b.block for i, b in enumerate(self._scheduled)}

    def get_scheduled_blocks(self) -> list[ScheduledScanBlock]:
        """Return the list of scheduled blocks with full timing metadata.

        Returns
        -------
        list[ScheduledScanBlock]
            All scheduled blocks in chronological order.
        """
        return getattr(self, '_scheduled', [])

    def print_schedule(self):
        """Print the schedule in a human-readable format to stdout."""
        if not (s := self.get_scheduled_blocks()):
            return print("No schedule. Run schedule() first.")
        print("\n" + "="*70 + "\nVLBI OBSERVATION SCHEDULE\n" + "="*70)
        for i, b in enumerate(s):
            print(f"\n{i+1}. {b.name}\n   {b.start_time.iso} - {b.end_time.iso}\n   "
                  f"{b.duration.value:.1f}min | {b.n_antennas} ant | {b.mean_elevation:.1f}° | {len(b.scans)} scans")
        print("\n" + "="*70)
