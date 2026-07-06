from typing import Optional
import numpy as np
from astropy import units as u
from astropy.time import Time
from astroplan import Observer, FixedTarget
from astroplan import Constraint, AltitudeConstraint, SunSeparationConstraint, \
                      MoonSeparationConstraint, min_best_rescale

"""Module that defines VLBI-specific `astroplan.Constraint` subclasses (elevation, azimuth,
hour angle, declination, and a per-station azimuth-dependent horizon), plus the
`horizon_min_elevation` helper used to evaluate the horizon constraint.
"""

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

    Parameters
    ----------
    azimuth_deg : np.ndarray
        Azimuth(s) of the target in degrees (any shape).
    horizon_az_deg : np.ndarray
        Azimuth sample points of the horizon, in degrees, sorted ascending within [0, 360].
    horizon_el_deg : np.ndarray
        Elevation (horizon) value in degrees at each azimuth sample point.

    Returns
    -------
    np.ndarray
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
        """Initializes the horizon constraint from its (azimuth, elevation) sample points.

        Parameters
        ----------
        horizon_az : `~astropy.units.Quantity`
            Azimuth sample points of the horizon (sorted ascending within [0, 360] degrees).
        horizon_el : `~astropy.units.Quantity`
            Horizon (minimum) elevation at each azimuth sample point.
        boolean_constraint : bool, optional
            If True (default), `compute_constraint` returns a boolean array (observable or not).
            If False, it returns a rescaled score usable for optimization instead of a hard cut.
        """
        self.horizon_az_deg: np.ndarray = np.atleast_1d(horizon_az.to(u.deg).value)
        self.horizon_el_deg: np.ndarray = np.atleast_1d(horizon_el.to(u.deg).value)
        self.boolean_constraint: bool = boolean_constraint

    def compute_constraint(self, times: Time, observer: Observer, targets: FixedTarget):
        """Evaluates the horizon constraint for the given times, observer, and targets.

        Parameters
        ----------
        times : `~astropy.time.Time`
            Times at which to evaluate the constraint.
        observer : `~astroplan.Observer`
            Observer (station) location used to compute the target's alt/az.
        targets : `~astroplan.FixedTarget`
            Target(s) to evaluate.

        Returns
        -------
        np.ndarray
            If `boolean_constraint` is True, a boolean array marking where the target elevation
            is above the local horizon. Otherwise, a non-negative array with the elevation margin
            above the horizon (0 where the target is below the horizon).
        """
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
    """Constraint on the elevation of the source. A thin alias of `astroplan.AltitudeConstraint`
    with VLBI-friendly naming; behavior (min/max, boolean_constraint) is inherited unchanged.
    """
    pass


class AzimuthConstraint(Constraint):
    """Constraint the azimuth for the station.
    """

    def __init__(self, min: Optional[u.Quantity] = None, max: Optional[u.Quantity] = None,
                 boolean_constraint: bool = True):
        """Initializes the azimuth constraint.

        Parameters
        ----------
        min : `~astropy.units.Quantity` or `None`, optional
            Minimum acceptable azimuth for the target. `None` indicates no limit (defaults to -10 deg).
        max : `~astropy.units.Quantity` or `None`, optional
            Maximum acceptable azimuth for the target. `None` indicates no limit (defaults to 370 deg).
        boolean_constraint : bool, optional
            If True (default), `compute_constraint` returns a boolean array (observable or not).
            If False, it returns a rescaled score usable for optimization instead of a hard cut.
        """
        self.min: u.Quantity = min if min is not None else -10*u.deg
        self.max: u.Quantity = max if max is not None else 370*u.deg
        self.boolean_constraint: bool = boolean_constraint

    def compute_constraint(self, times: Time, observer: Observer, targets: FixedTarget):
        """Evaluates the azimuth constraint for the given times, observer, and targets.

        Parameters
        ----------
        times : `~astropy.time.Time`
            Times at which to evaluate the constraint.
        observer : `~astroplan.Observer`
            Observer (station) location used to compute the target's azimuth.
        targets : `~astroplan.FixedTarget`
            Target(s) to evaluate.

        Returns
        -------
        np.ndarray
            If `boolean_constraint` is True, a boolean array marking where the azimuth falls
            within [min, max]. Otherwise, a rescaled score from `min_best_rescale`.
        """
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
        """Initializes the hour angle constraint.

        Parameters
        ----------
        min : `~astropy.units.Quantity` or `None`, optional
            Minimum acceptable hour angle for the target. `None` indicates no limit
            (defaults to -370 deg).
        max : `~astropy.units.Quantity` or `None`, optional
            Maximum acceptable hour angle for the target. `None` indicates no limit
            (defaults to 370 deg).
        boolean_constraint : bool, optional
            If True (default), `compute_constraint` returns a boolean array (observable or not).
            If False, it returns a rescaled score usable for optimization instead of a hard cut.
        """
        self.min = min if min is not None else -370*u.deg
        self.max = max if max is not None else 370*u.deg
        self.boolean_constraint = boolean_constraint

    def compute_constraint(self, times: Time, observer: Observer, targets: FixedTarget):
        """Evaluates the hour angle constraint for the given times, observer, and targets.

        The comparison is also repeated with min/max shifted by +24 hourangle, so that a wrapped
        range straddling 0/24h hour angle is still handled correctly.

        Parameters
        ----------
        times : `~astropy.time.Time`
            Times at which to evaluate the constraint.
        observer : `~astroplan.Observer`
            Observer (station) location used to compute the target's hour angle.
        targets : `~astroplan.FixedTarget`
            Target(s) to evaluate.

        Returns
        -------
        np.ndarray
            If `boolean_constraint` is True, a boolean array marking where the hour angle falls
            within [min, max] (accounting for the 24h wrap). Otherwise, a rescaled score from
            `min_best_rescale`.
        """
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
        """Initializes the declination constraint.

        Parameters
        ----------
        min : `~astropy.units.Quantity` or `None`, optional
            Minimum acceptable declination for the target. `None` indicates no limit
            (defaults to -90 deg).
        max : `~astropy.units.Quantity` or `None`, optional
            Maximum acceptable declination for the target. `None` indicates no limit
            (defaults to 91 deg).
        boolean_constraint : bool, optional
            If True (default), `compute_constraint` returns a boolean array (observable or not).
            If False, it returns a rescaled score usable for optimization instead of a hard cut.
        """
        self.min = min if min is not None else -90*u.deg
        self.max = max if max is not None else 91*u.deg
        self.boolean_constraint = boolean_constraint

    def compute_constraint(self, times: Time, observer: Observer, targets: FixedTarget):
        """Evaluates the declination constraint for the given targets.

        Note that `times` and `observer` are accepted for interface compatibility with
        `astroplan.Constraint` but are not used, since declination does not depend on time
        or observer location.

        Parameters
        ----------
        times : `~astropy.time.Time`
            Unused. Present for signature compatibility with `astroplan.Constraint`.
        observer : `~astroplan.Observer`
            Unused. Present for signature compatibility with `astroplan.Constraint`.
        targets : `~astroplan.FixedTarget`
            Target(s) to evaluate.

        Returns
        -------
        np.ndarray
            If `boolean_constraint` is True, a boolean array marking where the declination falls
            within [min, max]. Otherwise, a rescaled score from `min_best_rescale`.
        """
        target_declination = targets.dec
        if self.boolean_constraint:
            return ((self.min < target_declination) & (target_declination < self.max))
        else:
            rescale = min_best_rescale(target_declination, self.min, self.max, less_than_min=0)
            return rescale
