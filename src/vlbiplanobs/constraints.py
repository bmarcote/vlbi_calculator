from typing import Optional
import numpy as np
from astropy import units as u
from astropy.time import Time
from astroplan import Observer, FixedTarget
from astroplan import Constraint, AltitudeConstraint, SunSeparationConstraint, \
                      MoonSeparationConstraint, min_best_rescale

__all__: list[str] = ['ElevationConstraint', 'SunSeparationConstraint', 'MoonSeparationConstraint',
                      'AzimuthConstraint', 'HourAngleConstraint', 'DeclinationConstraint',
                      'HorizonConstraint', 'horizon_min_elevation']


def horizon_min_elevation(azimuth_deg: np.ndarray, horizon_az_deg: np.ndarray,
                          horizon_el_deg: np.ndarray) -> np.ndarray:
    """Returns the minimum observable elevation (the local horizon) at each given azimuth.

    The horizon is defined by a set of (azimuth, elevation) sample points. The elevation limit
    at an arbitrary azimuth is obtained by linear interpolation between the sample points, with
    the azimuth axis wrapping around at 360 degrees. A single sample point defines a flat horizon
    (constant elevation limit at all azimuths).

    Inputs
    - azimuth_deg : np.ndarray
        Azimuth(s) of the target in degrees (any shape).
    - horizon_az_deg : np.ndarray
        Azimuth sample points of the horizon, in degrees, sorted ascending within [0, 360].
    - horizon_el_deg : np.ndarray
        Elevation (horizon) value in degrees at each azimuth sample point.

    Returns
    - np.ndarray
        Minimum observable elevation in degrees at each input azimuth, same shape as `azimuth_deg`.
    """
    az = np.mod(np.asarray(azimuth_deg, dtype=np.float64), 360.0)
    xp = np.asarray(horizon_az_deg, dtype=np.float64)
    fp = np.asarray(horizon_el_deg, dtype=np.float64)
    if xp.size == 1:
        return np.full_like(az, fp[0])
    return np.interp(az, xp, fp, period=360.0)


class HorizonConstraint(Constraint):
    """Constraint enforcing an azimuth-dependent local horizon (minimum elevation).

    Some antennas are blocked by terrain or structures down to elevations higher than their
    nominal axis limit. This constraint models that local horizon as a set of (azimuth, elevation)
    sample points and requires the target elevation to be above the interpolated horizon at the
    target azimuth.
    """

    def __init__(self, horizon_az: u.Quantity, horizon_el: u.Quantity,
                 boolean_constraint: bool = True):
        """
        horizon_az : `~astropy.units.Quantity`
            Azimuth sample points of the horizon (sorted ascending within [0, 360] degrees).
        horizon_el : `~astropy.units.Quantity`
            Horizon (minimum) elevation at each azimuth sample point.
        """
        self.horizon_az_deg: np.ndarray = np.atleast_1d(horizon_az.to(u.deg).value)
        self.horizon_el_deg: np.ndarray = np.atleast_1d(horizon_el.to(u.deg).value)
        self.boolean_constraint: bool = boolean_constraint

    def compute_constraint(self, times: Time, observer: Observer, targets: FixedTarget):
        altaz = observer.altaz(times, targets)
        elevations = altaz.alt.to(u.deg).value
        azimuths = altaz.az.to(u.deg).value
        min_el = horizon_min_elevation(azimuths, self.horizon_az_deg, self.horizon_el_deg)
        if self.boolean_constraint:
            return elevations > min_el
        else:
            return np.clip(elevations - min_el, 0.0, None)

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
        self.min: u.Quantity = min if min is not None else -10*u.deg
        self.max: u.Quantity = max if max is not None else 370*u.deg
        self.boolean_constraint: bool = boolean_constraint

    def compute_constraint(self, times: Time, observer: Observer, targets: FixedTarget):
        azimuths = observer.altaz(times, targets).az
        if self.boolean_constraint:
            return ((self.min < azimuths) & (azimuths < self.max))
        else:
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
            return ((self.min < hour_angles) & (hour_angles < self.max)) | \
                   (((self.min + 24*u.hourangle) < hour_angles) & (hour_angles < (self.max + 24.0*u.hourangle)))
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
        target_declination = targets.dec
        if self.boolean_constraint:
            return ((self.min < target_declination) & (target_declination < self.max))
        else:
            rescale = min_best_rescale(target_declination, self.min, self.max, less_than_min=0)
            return rescale
