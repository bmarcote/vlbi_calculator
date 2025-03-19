import asyncio
from typing import Optional, Union, Iterable, Sequence, Self
from importlib import resources
from functools import partial
from itertools import product
import subprocess
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import numpy as np
import yaml
from rich import print as rprint
# import scipy.ndimage
import datetime as dt
from enum import Enum, auto
from dataclasses import dataclass
from astropy import units as u
from astropy import coordinates as coord
from astropy.time import Time, TimeDelta
from ortools.sat.python import cp_model
from .stations import Station, Stations
from .sources import Source, Scan, ScanBlock, SourceType, SourceNotVisible
from .observation import Observation


def _t2mjd10(time: Time) -> int:
    """Converts an astropy.time.Time object into a decimal MJD integer value, meaning int(MJD*10).
    """
    return int(time.mjd * 10)


def _t2time(time: int) -> Time:
    """Converts back a int(MJD*10) time stampt into a astropy.time.Time object.
    """
    return Time(time/10.0, format='mjd').datetime


class ScanBlockScheduler:
    def __init__(self, observation: Observation):
        """
        Args:
            observation (Observation): The observation to be schedued.
        """
        self.observation = observation
        self.starttime = _t2mjd10(observation.times[0])
        self.endtime = _t2mjd10(observation.times[-1])
        self.max_duration = self.endtime - self.starttime
        self.model = cp_model.CpModel()
        self.solver = cp_model.CpSolver()
        self.solution = None

    def create_variables(self):
        self.start_times = {}
        self.durations = {}
        self.end_times = {}

        for i, block in self.observation.scans.items():
            self.start_times[i] = self.model.NewIntVar(self.starttime, self.endtime,
                                                       f'start_time_{i}')
            self.durations[i] = self.model.NewIntVar(int(10*np.sum([s.duration.to(u.h).value for s in block.scans])),
                                                     self.max_duration, f'duration_{i}')
            self.end_times[i] = self.model.NewIntVar(self.starttime, self.endtime,
                                                     f'end_time_{i}')

            # Link start, duration, and end times
            self.model.Add(self.end_times[i] == self.start_times[i] + self.durations[i])

    def add_constraints(self):
        # Ensure all scan blocks are within observation window
        for i in self.observation.scans.keys():
            self.model.Add(self.start_times[i] >= self.starttime)
            self.model.Add(self.end_times[i] <= self.endtime)

        # Non-overlapping constraint
        for ii, i in enumerate(self.observation.scans.keys()):
            for j in list(self.observation.scans.keys())[ii+1:]:
                self.model.AddNoOverlap([
                    self.model.NewIntervalVar(self.start_times[i], self.durations[i], self.end_times[i],
                                              f'interval_{i}'),
                    self.model.NewIntervalVar(self.start_times[j], self.durations[j], self.end_times[j],
                                              f'interval_{j}')
                ])

    def optimize_station_participation(self):
        # BUG: I don't think this works as it is!
        for i, block in self.observation.scans.items():
            participation_var = self.model.NewIntVar(0, len(self.stations), f'participation_{i}')
            self.model.Add(participation_var == sum(
                           self.observation.is_observable_at(self.self.start_times[i]).values()
                           ))
            self.model.Maximize(participation_var)

    def optimize_coordinate_separation(self):
        # BUG: this also doesn't work.  How I compute one block and the next one in this context?
        for i, block in self.observation.scans.items():
            separation = self.model.NewIntVar(0, 360, f'separation_{i}')
            self.model.Add(separation == int(abs(
                self.observation.scans[i].sources(SourceType.TARGET)[0].coordinates -
                self.observation.scans[i+1].sources(SourceType.TARGET)[0].coordinates
            )))
            self.model.Minimize(separation)

    def optimize_elevation(self):
        time_index = lambda i: np.where(self.observation.times.mjd == min(self.observation.times.mjd,
                                        key=lambda x: abs(x-self.start_times[i]/10)))
        for i, elevations in enumerate(self.observation.elevations()):
            for station, elev in self.observation.elevations()[elevations].items():
                elevation = self.model.NewIntVar(0, 90, f'elevation_{i}_{station}')
                self.model.Add(elevation == elev[time_index(i)])
                self.model.Add(elevation >= 20)
                self.model.Add(elevation <= 82)

    def solve(self):
        self.create_variables()
        self.add_constraints()
        # self.optimize_station_participation()
        # self.optimize_elevation()
        # self.optimize_coordinate_separation()

        status = self.solver.Solve(self.model)

        if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
            schedule = []
            for i, block in self.observation.scans.items():
                start = self.solver.Value(self.start_times[i])
                duration = self.solver.Value(self.durations[i])
                schedule.append((i, start, duration))

            self.solution = schedule
        else:
            self.solution = None

    def print_schedule(self):
        for block, start, duration in self.solution:
            print(f"ScanBlock {block}: Start={_t2time(start)}, Duration={duration}")


# scheduler = ScanBlockScheduler(scan_blocks, stations, observation_start, observation_end)
# optimal_schedule = scheduler.solve()
#
# if optimal_schedule:
#     for block, start, duration in optimal_schedule:
#         print(f"ScanBlock {block.id}: Start={start}, Duration={duration}")
# else:
#     print("No feasible schedule found.")
