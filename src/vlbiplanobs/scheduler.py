# -*- coding: utf-8 -*-
# Licensed under GPLv3+ - see LICENSE
"""VLBI Observation Scheduler using Constraint Programming.

This module provides a scheduler that uses OR-Tools CP-SAT solver to find optimal
schedules for VLBI observations with multiple scan blocks. The scheduler optimizes for:
- Maximum antenna participation (more antennas observing each block)
- Maximum elevation (better signal quality)
- Minimum slewing time between blocks (minimize time lost to antenna movement)
"""

from dataclasses import dataclass
from typing import Optional
import numpy as np
from astropy import units as u
from astropy.time import Time
from ortools.sat.python import cp_model
from .sources import SourceType
from .observation import Observation


@dataclass
class ScheduledBlock:
    """Represents a scheduled scan block with timing information.

    Attributes
    - name : str
        Name/identifier of the scan block.
    - start_time : Time
        Start time of the scheduled block.
    - end_time : Time
        End time of the scheduled block.
    - duration_minutes : float
        Duration of the block in minutes.
    - n_antennas : int
        Number of antennas that can observe during this block.
    - mean_elevation : float
        Mean elevation across all observing antennas (degrees).
    """
    name: str
    start_time: Time
    end_time: Time
    duration_minutes: float
    n_antennas: int
    mean_elevation: float


class ScanBlockScheduler:
    """Scheduler for VLBI scan blocks using constraint programming.

    Uses OR-Tools CP-SAT solver to find optimal schedules that maximize antenna
    participation and elevation while minimizing slewing time between blocks.
    """

    def __init__(self, observation: Observation, min_antennas: int = 2,
                 min_elevation: float = 10.0, max_elevation: float = 85.0):
        """Initialize the scheduler with an observation.

        Inputs
        - observation : Observation
            The VLBI observation containing scan blocks to schedule.
        - min_antennas : int
            Minimum number of antennas required to observe each block (default: 2).
        - min_elevation : float
            Minimum elevation angle in degrees for valid observations (default: 10.0).
        - max_elevation : float
            Maximum elevation angle in degrees for valid observations (default: 85.0).
        """
        self.observation = observation
        self.min_antennas = min_antennas
        self.min_elevation = min_elevation
        self.max_elevation = max_elevation
        self.model = cp_model.CpModel()
        self.solver = cp_model.CpSolver()
        self.solution: Optional[list[ScheduledBlock]] = None
        self._block_names = list(observation.scans.keys())
        self._n_blocks = len(self._block_names)
        self._n_times = len(observation.times)
        self._n_antennas = len(observation.stations)
        self._precompute_data()

    def _precompute_data(self):
        """Pre-compute visibility and elevation data for use in CP model.

        CP-SAT requires integer arithmetic, so we pre-compute lookup tables for
        visibility (boolean) and elevation (scaled to integers).
        """
        is_obs = self.observation.is_observable()
        self._visibility = np.zeros((self._n_blocks, self._n_times), dtype=int)
        for bi, block_name in enumerate(self._block_names):
            if block_name in is_obs:
                self._visibility[bi, :] = np.sum(np.array(list(is_obs[block_name].values())), axis=0)

        elevs = self.observation.elevations()
        self._elevation = np.zeros((self._n_blocks, self._n_times), dtype=int)
        for bi, block_name in enumerate(self._block_names):
            block = self.observation.scans[block_name]
            targets = block.sources(SourceType.TARGET)
            if not targets:
                targets = block.sources()
            if targets and targets[0].name in elevs:
                src_elevs = elevs[targets[0].name]
                mean_elev = np.nanmean(np.array([src_elevs[ant].value for ant in src_elevs]), axis=0)
                self._elevation[bi, :] = np.where(np.isnan(mean_elev), 0, (mean_elev * 10).astype(int))

        self._separation = np.zeros((self._n_blocks, self._n_blocks), dtype=int)
        for bi, name_i in enumerate(self._block_names):
            block_i = self.observation.scans[name_i]
            targets_i = block_i.sources(SourceType.TARGET) or block_i.sources()
            if not targets_i:
                continue
            coord_i = targets_i[0].coord
            for bj, name_j in enumerate(self._block_names):
                if bi == bj:
                    continue
                block_j = self.observation.scans[name_j]
                targets_j = block_j.sources(SourceType.TARGET) or block_j.sources()
                if not targets_j:
                    continue
                coord_j = targets_j[0].coord
                sep = coord_i.separation(coord_j).arcmin
                self._separation[bi, bj] = int(sep)

        time_step = (self.observation.times[1] - self.observation.times[0]).to(u.min).value
        self._time_step_minutes = time_step
        self._block_durations = {}
        for bi, block_name in enumerate(self._block_names):
            block = self.observation.scans[block_name]
            total_min = sum(s.duration.to(u.min).value for s in block.scans)
            self._block_durations[bi] = max(1, int(np.ceil(total_min / time_step)))

    def _create_variables(self):
        """Create CP-SAT variables for the scheduling problem."""
        self._start_vars = {}
        self._end_vars = {}
        self._interval_vars = {}
        self._order_vars = {}

        for bi in range(self._n_blocks):
            duration = self._block_durations[bi]
            max_start = self._n_times - duration

            self._start_vars[bi] = self.model.NewIntVar(0, max_start, f'start_{bi}')
            self._end_vars[bi] = self.model.NewIntVar(duration, self._n_times, f'end_{bi}')
            self.model.Add(self._end_vars[bi] == self._start_vars[bi] + duration)
            self._interval_vars[bi] = self.model.NewIntervalVar(
                self._start_vars[bi], duration, self._end_vars[bi], f'interval_{bi}')

        for bi in range(self._n_blocks):
            self._order_vars[bi] = self.model.NewIntVar(0, self._n_blocks - 1, f'order_{bi}')

        self._before_vars = {}
        for bi in range(self._n_blocks):
            for bj in range(self._n_blocks):
                if bi != bj:
                    self._before_vars[(bi, bj)] = self.model.NewBoolVar(f'before_{bi}_{bj}')

    def _add_constraints(self):
        """Add constraints to the CP model."""
        self.model.AddNoOverlap(list(self._interval_vars.values()))
        self.model.AddAllDifferent(list(self._order_vars.values()))

        for bi in range(self._n_blocks):
            for bj in range(self._n_blocks):
                if bi != bj:
                    self.model.Add(self._order_vars[bi] < self._order_vars[bj]).OnlyEnforceIf(
                        self._before_vars[(bi, bj)])
                    self.model.Add(self._order_vars[bi] >= self._order_vars[bj]).OnlyEnforceIf(
                        self._before_vars[(bi, bj)].Not())
                    self.model.Add(self._end_vars[bi] <= self._start_vars[bj]).OnlyEnforceIf(
                        self._before_vars[(bi, bj)])

        for bi in range(self._n_blocks):
            valid_starts = [ti for ti in range(self._n_times - self._block_durations[bi])
                            if self._visibility[bi, ti] >= self.min_antennas]
            if valid_starts:
                self.model.AddAllowedAssignments([self._start_vars[bi]], [[t] for t in valid_starts])
            else:
                print(f"Warning: Block {self._block_names[bi]} has no valid start times with "
                      f"{self.min_antennas}+ antennas. Relaxing constraint.")

    def _add_objectives(self, weight_antennas: int = 100, weight_elevation: int = 10,
                        weight_slewing: int = 1):
        """Add optimization objectives to the CP model.

        Inputs
        - weight_antennas : int
            Weight for antenna participation objective (higher = more important).
        - weight_elevation : int
            Weight for elevation objective (higher = more important).
        - weight_slewing : int
            Weight for slewing minimization (higher = more important).
        """
        objective_terms = []

        for bi in range(self._n_blocks):
            antenna_score = self.model.NewIntVar(0, self._n_antennas * 10, f'ant_score_{bi}')
            self.model.AddElement(self._start_vars[bi], self._visibility[bi].tolist(), antenna_score)
            objective_terms.append(weight_antennas * antenna_score)

        for bi in range(self._n_blocks):
            elev_score = self.model.NewIntVar(0, 900, f'elev_score_{bi}')
            self.model.AddElement(self._start_vars[bi], self._elevation[bi].tolist(), elev_score)
            objective_terms.append(weight_elevation * elev_score)

        for bi in range(self._n_blocks):
            for bj in range(self._n_blocks):
                if bi != bj:
                    separation_penalty = self.model.NewIntVar(0, 360 * 60, f'sep_penalty_{bi}_{bj}')
                    self.model.Add(separation_penalty == self._separation[bi, bj]).OnlyEnforceIf(
                        self._before_vars[(bi, bj)])
                    self.model.Add(separation_penalty == 0).OnlyEnforceIf(
                        self._before_vars[(bi, bj)].Not())
                    objective_terms.append(-weight_slewing * separation_penalty)

        self.model.Maximize(sum(objective_terms))

    def solve(self, time_limit_seconds: float = 60.0, weight_antennas: int = 100,
              weight_elevation: int = 10, weight_slewing: int = 1) -> Optional[list[ScheduledBlock]]:
        """Solve the scheduling problem and return the optimal schedule.

        Inputs
        - time_limit_seconds : float
            Maximum time to spend solving (default: 60 seconds).
        - weight_antennas : int
            Weight for antenna participation objective.
        - weight_elevation : int
            Weight for elevation objective.
        - weight_slewing : int
            Weight for slewing minimization objective.

        Returns
        - list[ScheduledBlock] or None
            List of scheduled blocks in temporal order, or None if no solution found.
        """
        self._create_variables()
        self._add_constraints()
        self._add_objectives(weight_antennas, weight_elevation, weight_slewing)

        self.solver.parameters.max_time_in_seconds = time_limit_seconds
        status = self.solver.Solve(self.model)

        if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            self.solution = self._extract_solution()
            return self.solution
        else:
            self.solution = None
            return None

    def _extract_solution(self) -> list[ScheduledBlock]:
        """Extract the solution from the solver and return ordered ScheduledBlock list."""
        blocks_with_order = []
        for bi in range(self._n_blocks):
            start_idx = self.solver.Value(self._start_vars[bi])
            end_idx = self.solver.Value(self._end_vars[bi])
            order = self.solver.Value(self._order_vars[bi])

            start_time = self.observation.times[start_idx]
            end_time = self.observation.times[min(end_idx, self._n_times - 1)]
            duration = self._block_durations[bi] * self._time_step_minutes
            n_antennas = int(self._visibility[bi, start_idx])
            mean_elev = self._elevation[bi, start_idx] / 10.0

            scheduled = ScheduledBlock(
                name=self._block_names[bi],
                start_time=start_time,
                end_time=end_time,
                duration_minutes=duration,
                n_antennas=n_antennas,
                mean_elevation=mean_elev
            )
            blocks_with_order.append((order, scheduled))

        blocks_with_order.sort(key=lambda x: x[0])
        return [block for _, block in blocks_with_order]

    def print_schedule(self):
        """Print the schedule in a human-readable format."""
        if self.solution is None:
            print("No schedule available. Run solve() first.")
            return

        print("\n" + "=" * 70)
        print("VLBI OBSERVATION SCHEDULE")
        print("=" * 70)
        total_slew = 0.0
        for i, block in enumerate(self.solution):
            print(f"\n{i + 1}. {block.name}")
            print(f"   Start: {block.start_time.iso}")
            print(f"   End:   {block.end_time.iso}")
            print(f"   Duration: {block.duration_minutes:.1f} min")
            print(f"   Antennas: {block.n_antennas}")
            print(f"   Mean elevation: {block.mean_elevation:.1f}°")
            if i > 0:
                prev_block = self.solution[i - 1]
                bi = self._block_names.index(prev_block.name)
                bj = self._block_names.index(block.name)
                slew = self._separation[bi, bj] / 60.0
                total_slew += slew
                print(f"   Slew from previous: {slew:.1f}°")

        print("\n" + "-" * 70)
        print(f"Total angular slewing: {total_slew:.1f}°")
        print("=" * 70)
