# -*- coding: utf-8 -*-
# Licensed under GPLv3+ - see LICENSE
# from __future__ import annotations
# from collections import abc
from typing import Optional
# import configparser
# from importlib import resources
# import numpy as np
from astropy import units as u
# from astropy import coordinates as coord
# # from astropy.io import ascii
from astropy.time import Time
from astroplan import Observer, FixedTarget
from astroplan import Constraint, AltitudeConstraint, SunSeparationConstraint, MoonSeparationConstraint, \
                      min_best_rescale
# from dataclasses import dataclass
# from enum import Enum, auto


# Limits for the station to check if it can observe, adapted for VLBI
# - SunSeparationConstraint(min, max=None)
# - AltitudeConstraint(min, max)
# - MoonSeparationConstraint(min, max=None)

class ElevationConstraint(AltitudeConstraint):
    """Constraint on the elevation of the source.
    """
    pass


class AzimuthConstraint(Constraint):
    """Constraint the azimuth for the station
    """
    def __init__(self, min: Optional[u.Quantity] = None, max: Optional[u.Quantity] = None,
                 boolean_constraint: bool = True):
        """
        min : `~astropy.units.Quantity` or `None` (optional)
            Minimum acceptable azimuth for the target. `None` indicates no limit.
        max : `~astropy.units.Quantity` or `None` (optional)
            Minimum acceptable azimuth for the target. `None` indicates no limit.
        """
        self.min = min if min is not None else -10*u.deg
        self.max = max if max is not None else 370*u.deg
        self.boolean_constraint = boolean_constraint

    def compute_constraint(self, times: Time, observer: Observer, targets: FixedTarget):
        azimuths = observer.altaz(times, targets).az
        if self.boolean_constraint:
            return ((self.min < azimuths) & (azimuths < self.max))
        # if we want to return a non-boolean score
        else:
            # rescale the vega_separation values so that they become
            # scores between zero and one
            rescale = min_best_rescale(azimuths, self.min, self.max, less_than_min=0)
            return rescale


class HourAngleConstraint(Constraint):
    """Constraint the hour angle for the station.
    """
    def __init__(self, min: Optional[u.Quantity] = None, max: Optional[u.Quantity] = None,
                 boolean_constraint: bool = True):
        """
        min : `~astropy.units.Quantity` or `None` (optional)
            Minimum acceptable hour angle for the target. `None` indicates no limit.
        max : `~astropy.units.Quantity` or `None` (optional)
            Maximum acceptable hour angle for the target. `None` indicates no limit.
        """
        self.min = min if min is not None else -370*u.deg
        self.max = max if max is not None else 370*u.deg
        self.boolean_constraint = boolean_constraint

    def compute_constraint(self, times: Time, observer: Observer, targets: FixedTarget):
        hour_angles = observer.target_hour_angle(times, targets)
        if self.boolean_constraint:
            return ((self.min < hour_angles) & (hour_angles < self.max))
        else:
            rescale = min_best_rescale(hour_angles, self.min, self.max, less_than_min=0)
            return rescale


class DeclinationConstraint(Constraint):
    """Constraint the source declination for the station.
    """
    def __init__(self, min: Optional[u.Quantity] = None, max: Optional[u.Quantity] = None,
                 boolean_constraint: bool = True):
        """
        min : `~astropy.units.Quantity` or `None` (optional)
            Minimum acceptable declination for the target. `None` indicates no limit.
        max : `~astropy.units.Quantity` or `None` (optional)
            Maximum acceptable declination for the target. `None` indicates no limit.
        """
        self.min = min if min is not None else -90*u.deg
        self.max = max if max is not None else 91*u.deg
        self.boolean_constraint = boolean_constraint

    def compute_constraint(self, times: Time, observer: Observer, targets: FixedTarget):
        target_declination = targets.coord.dec
        if self.boolean_constraint:
            return ((self.min < target_declination) & (target_declination < self.max))
        else:
            rescale = min_best_rescale(target_declination, self.min, self.max, less_than_min=0)
            return rescale




