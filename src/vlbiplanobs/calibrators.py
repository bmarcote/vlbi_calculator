"""Calibrator source module for finding fringe finders and nearby sources."""

import argparse
import json
import logging
import functools
from typing import NamedTuple, Optional, Self
from importlib import resources
import numpy as np
import erfa
from astropy import units as u, coordinates as coord
from astropy.time import Time

from rich import print as rprint, box
from rich.table import Table
from rich_argparse import RawTextRichHelpFormatter
from urllib import parse

from .sources import Source, SourceType, SourceCatalog
from .stations import Stations, MountType
from . import observation as obs

_log = logging.getLogger(__name__)

_RFC_BANDS: dict[str, str] = {'l': 's', 's': 's', 'c': 'c', 'm': 'c', 'x': 'x', 'u': 'u', 'k': 'k', 'q': 'k'}
_DEFAULT_MIN_ELEVATION: u.Quantity = 20 * u.deg
_DEFAULT_MIN_FLUX: u.Quantity = 1.0 * u.Jy
_BAND_INDEX: dict[str, int] = {'s': 0, 'c': 1, 'x': 2, 'u': 3, 'k': 4}
_WAVELENGTH_BANDS: dict[str, str] = {'18cm': 's', '21cm': 's', '13cm': 'c', '6cm': 'c', '5cm': 'c',
                       '3.6cm': 'x', '2cm': 'u', '1.3cm': 'k', '0.7cm': 'k'}


def _round_to_nearest_wavelength(band_str: str) -> str:
    """Round an unknown wavelength to the nearest known wavelength.

    Parameters
    ----------
    band_str : str
        The wavelength string to round.

    Returns
    -------
    str
        The nearest known wavelength string.
    """
    if band_str in _WAVELENGTH_BANDS:
        return band_str
    if not band_str.lower().endswith('cm'):
        return band_str
    return min(_WAVELENGTH_BANDS.keys(), key=lambda x: abs(float(x[:-2]) - float(band_str[:-2])))


class CalibratorSource(Source):
    """Represents a calibrator source from the RFC catalog.

    Attributes
    ----------
    ivsname : str
        IVS name of the source.
    n_observations : int
        Number of observations.
    flux_resolved : np.ndarray
        Resolved flux values per band.
    flux_unresolved : np.ndarray
        Unresolved flux values per band.
    is_calibrator : bool
        Whether the source is a calibrator.
    """
    __slots__ = ('ivsname', 'n_observations', 'flux_resolved', 'flux_unresolved', 'is_calibrator',
                 '_ra_deg', '_dec_deg', '_coord')

    def __init__(self, name: str, ivsname: str, ra_deg: float, dec_deg: float, n_observations: int,
                 flux_resolved: np.ndarray, flux_unresolved: np.ndarray, is_calibrator: bool):
        """Initializes a calibrator source from already-parsed RFC catalog values.

        The astropy SkyCoord (`coord`) is built lazily on first access, because building ~13k SkyCoord
        objects eagerly dominated the RFCCatalog load time. Source.__init__ is intentionally not called;
        the attributes it would set are assigned explicitly here.

        Parameters
        ----------
        name : str
            J2000 name of the source.
        ivsname : str
            IVS name of the source.
        ra_deg, dec_deg : float
            ICRS right ascension and declination in degrees.
        n_observations : int
            Number of observations in the RFC catalog.
        flux_resolved, flux_unresolved : np.ndarray
            Flux values (Jy) per RFC band (s, c, x, u, k).
        is_calibrator : bool
            Whether the RFC flags the source as calibrator ('C').
        """
        self.name = name
        self._type = SourceType.PHASECAL
        self._flux = None
        self._notes = None
        self._other_names = []
        self._ra_deg: float = float(ra_deg)
        self._dec_deg: float = float(dec_deg)
        self._coord: Optional[coord.SkyCoord] = None
        self.ivsname: str = ivsname
        self.n_observations: int = n_observations
        self.flux_resolved: np.ndarray = flux_resolved
        self.flux_unresolved: np.ndarray = flux_unresolved
        self.is_calibrator: bool = is_calibrator

    @property
    def coord(self) -> "coord.SkyCoord":
        """ICRS coordinates of the source (built lazily from ra_deg/dec_deg)."""
        if self._coord is None:
            self._coord = coord.SkyCoord(ra=self._ra_deg * u.deg, dec=self._dec_deg * u.deg)
        return self._coord

    @coord.setter
    def coord(self, value: "coord.SkyCoord"):
        """Overrides the coordinates of the source (keeps ra_deg/dec_deg consistent)."""
        self._coord = value
        self._ra_deg = float(value.ra.deg)
        self._dec_deg = float(value.dec.deg)

    @property
    def ra_deg(self) -> float:
        """Right ascension in degrees."""
        return self._ra_deg

    @property
    def dec_deg(self) -> float:
        """Declination in degrees."""
        return self._dec_deg

    def unresolved_flux(self, band: str) -> float:
        """Get unresolved flux for a specific band.

        Parameters
        ----------
        band : str
            RFC band code (s/c/x/u/k).

        Returns
        -------
        float
            Unresolved flux in Jy, or 0.0 if band not available.
        """
        idx = _BAND_INDEX.get(band)
        return float(self.flux_unresolved[idx]) if idx is not None else 0.0

    def resolved_flux(self, band: str) -> float:
        """Get resolved flux for a specific band.

        Parameters
        ----------
        band : str
            RFC band code (s/c/x/u/k).

        Returns
        -------
        float
            Resolved flux in Jy, or 0.0 if band not available.
        """
        idx = _BAND_INDEX.get(band)
        return float(self.flux_resolved[idx]) if idx is not None else 0.0

    def get_flux_at_band(self, band: str) -> tuple[float, float]:
        """Get resolved and unresolved flux at a specific band with interpolation.

        Parameters
        ----------
        band : str
            RFC band code or wavelength string.

        Returns
        -------
        tuple[float, float]
            (resolved_flux, unresolved_flux) in Jy.
        """
        if band in _WAVELENGTH_BANDS:
            band = _WAVELENGTH_BANDS[band]
        elif 'cm' in band.lower():
            rounded_band = _round_to_nearest_wavelength(band)
            if rounded_band in _WAVELENGTH_BANDS:
                band = _WAVELENGTH_BANDS[rounded_band]

        idx = _BAND_INDEX.get(band)
        if idx is None:
            return 0.0, 0.0

        if idx < len(self.flux_unresolved) and self.flux_unresolved[idx] > 0:
            return float(self.flux_resolved[idx]), float(self.flux_unresolved[idx])

        return self._interpolate_flux(band)

    def _interpolate_flux(self, target_band: str) -> tuple[float, float]:
        """Interpolate flux for a band using nearby bands.

        Parameters
        ----------
        target_band : str
            Target RFC band code.

        Returns
        -------
        tuple[float, float]
            (resolved_flux, unresolved_flux) in Jy.
        """
        # Find the wavelength for the target band - check both direct band and wavelength mappings
        target_wavelength = None
        for wl, band in _WAVELENGTH_BANDS.items():
            if band == target_band:
                target_wavelength = float(wl[:-2])
                break

        # If no direct mapping found, try to find a representative wavelength for the band
        if target_wavelength is None and target_band in _BAND_INDEX:
            # Use a representative wavelength for each band
            representative_wavelengths = {'s': 18.0, 'c': 6.0, 'x': 3.6, 'u': 2.0, 'k': 1.3}
            target_wavelength = representative_wavelengths.get(target_band)

        if target_wavelength is None:
            return 0.0, 0.0

        # Vectorized approach: collect all available data at once
        available_data = []
        for band, idx in _BAND_INDEX.items():
            if idx < len(self.flux_unresolved) and self.flux_unresolved[idx] > 0:
                for wl, wl_band in _WAVELENGTH_BANDS.items():
                    if wl_band == band:
                        available_data.append((float(wl[:-2]), float(self.flux_resolved[idx]),
                                               float(self.flux_unresolved[idx])))
                        break

        if not available_data:
            return 0.0, 0.0

        if len(available_data) == 1:
            return available_data[0][1], available_data[0][2]

        wavelengths, resolved_fluxes, unresolved_fluxes = map(np.array, zip(*available_data))
        distances = np.abs(wavelengths - target_wavelength)
        closest_indices = np.argsort(distances)[:2]

        wl1, wl2 = wavelengths[closest_indices[0]], wavelengths[closest_indices[1]]
        res1, res2 = resolved_fluxes[closest_indices[0]], resolved_fluxes[closest_indices[1]]
        unres1, unres2 = unresolved_fluxes[closest_indices[0]], unresolved_fluxes[closest_indices[1]]

        if wl1 != wl2 and target_wavelength > 0:
            log_wl1, log_wl2 = np.log10(wl1), np.log10(wl2)
            log_target_wl = np.log10(target_wavelength)

            if res1 > 0 and res2 > 0:
                log_res1, log_res2 = np.log10(res1), np.log10(res2)
                log_res_interp = log_res1 + (log_res2 - log_res1) * (log_target_wl - \
                    log_wl1) / (log_wl2 - log_wl1)
                res_interp = 10**log_res_interp
            else:
                res_interp = (res1 + res2) / 2

            if unres1 > 0 and unres2 > 0:
                log_unres1, log_unres2 = np.log10(unres1), np.log10(unres2)
                log_unres_interp = log_unres1 + (log_unres2 - log_unres1) * (log_target_wl - log_wl1) / \
                    (log_wl2 - log_wl1)
                unres_interp = 10**log_unres_interp
            else:
                unres_interp = (unres1 + unres2) / 2
        else:
            res_interp = (res1 + res2) / 2
            unres_interp = (unres1 + unres2) / 2

        return res_interp, unres_interp

    def get_skycoord(self) -> "coord.SkyCoord":
        """Get the SkyCoord object for this source.

        Returns
        -------
        astropy.coordinates.SkyCoord
            The source coordinates.
        """
        return self.coord

    def get_astrogeo_link(self) -> str:
        """Generate astrogeo.org calibrator search link for this source.

        Returns
        -------
        str
            URL to astrogeo calibrator search.
        """
        # signed_dms keeps the sign separate, so -1 deg < Dec < 0 deg is not lost as '-0' degrees.
        ra_h, ra_m, ra_s = self.coord.ra.hms
        dec_sign, dec_d, dec_m, dec_s = self.coord.dec.signed_dms
        sign_char = '-' if dec_sign < 0 else '+'
        fmt = "ra={:02.0f}:{:02.0f}:{:06.3f}&dec={}{:02.0f}:{:02.0f}:{:06.3f}&num_sou=1&format=html"
        source_coord_str = parse.quote(fmt.format(ra_h, ra_m, ra_s, sign_char, dec_d, dec_m, dec_s), safe='=&')
        return f"https://astrogeo.org/cgi-bin/calib_search_form.csh?{source_coord_str}"

    def get_observed_bands(self) -> str:
        """Get comma-separated list of bands with observed flux.

        Returns
        -------
        str
            Comma-separated uppercase band letters (e.g., 'S,C,X').
        """
        return ','.join(band.upper()
                        for band, idx in _BAND_INDEX.items()
                        if self.flux_unresolved[idx] > 0 or self.flux_resolved[idx] > 0)


class _RFCRows(NamedTuple):
    """Immutable, parsed content of an RFC catalog file (one entry per valid source row).

    All numpy arrays are marked read-only, so they can be safely shared between RFCCatalog instances.
    """
    names: tuple[str, ...]
    ivsnames: tuple[str, ...]
    ra_deg: np.ndarray
    dec_deg: np.ndarray
    n_obs: np.ndarray
    flux: np.ndarray  # shape (n_sources, 5 bands, 2) with [:, :, 0]=resolved, [:, :, 1]=unresolved (float32, Jy)
    is_calibrator: np.ndarray


def _parse_rfc_line(line: str) -> Optional[tuple]:
    """Parses one RFC catalog row.

    Parameters
    ----------
    line : str
        Raw text line from the RFC catalog file.

    Returns
    -------
    tuple or None
        (name, ivsname, ra_deg, dec_deg, n_obs, flux_values(10 floats), is_calibrator),
        or None if the line is a comment/unreliable entry or cannot be parsed.
    """
    if not line or line[0] in '#U' or len(line) < 100:
        return None

    cols = line.split()
    if len(cols) < 25:
        return None

    try:
        flux_values = []
        for f_res, f_unres in zip(cols[13:23:2], cols[14:24:2]):
            flux_values.append(0.0 if f_res[0] == '<' else float(f_res))
            flux_values.append(0.0 if f_unres[0] == '<' else float(f_unres))

        ra_deg = (float(cols[3]) + float(cols[4]) / 60.0 + float(cols[5]) / 3600.0) * 15.0
        dec_deg = (1.0 if cols[6][0] != '-' else -1.0) * \
            (abs(float(cols[6])) + float(cols[7]) / 60.0 + float(cols[8]) / 3600.0)
    except (ValueError, IndexError):
        return None

    n_obs = int(cols[12]) if cols[12].lstrip('-').isdigit() else 0
    return cols[2], cols[1], ra_deg, dec_deg, n_obs, flux_values, cols[0] == 'C'


@functools.lru_cache(maxsize=8)
def _parse_rfc_file(catalog_path: str) -> _RFCRows:
    """Parses an RFC catalog file once and caches the immutable result (keyed by path).

    Parameters
    ----------
    catalog_path : str
        Path to the RFC catalog text file.

    Returns
    -------
    _RFCRows
        Parsed rows with read-only numpy arrays.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    RuntimeError
        If the file cannot be read/parsed.
    """
    try:
        with open(catalog_path, 'rt') as fin:
            parsed = [row for row in (_parse_rfc_line(line) for line in fin) if row is not None]
    except FileNotFoundError as e:
        raise FileNotFoundError(f'RFC catalog file not found: {catalog_path}') from e
    except (OSError, UnicodeDecodeError) as e:
        raise RuntimeError(f'Error reading RFC catalog file {catalog_path}: {e}') from e

    n_rows = len(parsed)
    rows = _RFCRows(names=tuple(r[0] for r in parsed), ivsnames=tuple(r[1] for r in parsed),
                    ra_deg=np.array([r[2] for r in parsed], dtype=np.float64),
                    dec_deg=np.array([r[3] for r in parsed], dtype=np.float64),
                    n_obs=np.array([r[4] for r in parsed], dtype=np.int64),
                    flux=np.array([r[5] for r in parsed], dtype=np.float32).reshape(n_rows, 5, 2),
                    is_calibrator=np.array([r[6] for r in parsed], dtype=bool))
    for arr in (rows.ra_deg, rows.dec_deg, rows.n_obs, rows.flux, rows.is_calibrator):
        arr.flags.writeable = False

    _log.info("Parsed RFC catalog %s: %d sources", catalog_path, n_rows)
    return rows


class RFCCatalog:
    """RFC (Radio Fundamental Catalog) of VLBI calibrator sources.

    The catalog file is parsed only once per path (module-level cache); each instance applies its own
    min_flux/band/include_missing filter and owns its own list of CalibratorSource objects.

    Attributes
    ----------
    _sources : list[CalibratorSource]
        List of calibrator sources.
    _min_flux : float
        Minimum flux threshold in Jy.
    _band : str
        Default RFC band code.
    _catalog_filename : str or None
        Path to catalog file, or None for default.
    _name_index : dict[str, CalibratorSource]
        Index mapping upper-case source (J2000) names to CalibratorSource objects.
    _ivsname_index : dict[str, CalibratorSource]
        Index mapping upper-case IVS names to CalibratorSource objects.
    _ra_arr : np.ndarray
        Array of source right ascensions in degrees.
    _dec_arr : np.ndarray
        Array of source declinations in degrees.
    _include_missing : bool
        Whether to include sources with missing band measurements.
    """
    __slots__ = ('_sources', '_min_flux', '_band', '_catalog_filename',
                 '_name_index', '_ivsname_index', '_ra_arr', '_dec_arr', '_include_missing')

    def __init__(self, catalog_filename: Optional[str] = None, min_flux: u.Quantity = _DEFAULT_MIN_FLUX,
                 band: str = 'c', include_missing: bool = False):
        """Loads the RFC catalog, keeping only the sources passing the flux filter.

        Parameters
        ----------
        catalog_filename : str, optional
            Path to the RFC catalog file. If None, the newest packaged RFC catalog is used.
        min_flux : astropy.units.Quantity or float, optional
            Minimum unresolved flux at `band` (Jy if a float). Default 1 Jy.
        band : str, optional
            RFC band code (s/c/x/u/k) used for the flux filter. Default 'c'.
        include_missing : bool, optional
            If True, keep sources without a measurement at `band` (phase-cal mode). Default False.

        Raises
        ------
        ValueError
            If `band` is not a valid RFC band code.
        FileNotFoundError
            If the catalog file is not found.
        RuntimeError
            If the catalog file cannot be read.
        """
        if band not in _BAND_INDEX:
            raise ValueError(f"Unknown RFC band '{band}'. Valid bands: {list(_BAND_INDEX)}.")

        self._min_flux = min_flux.to(u.Jy).value if hasattr(min_flux, 'to') else min_flux
        self._band = band
        self._include_missing = include_missing
        self._catalog_filename = catalog_filename
        self._load_catalog()

    def _get_catalog_path(self) -> str:
        """Get the path to the RFC catalog file.

        Returns
        -------
        str
            Path to the catalog file.
        """
        if self._catalog_filename is not None:
            return str(self._catalog_filename)
        rfc_files = tuple(r.name for r in resources.files('vlbiplanobs.data').iterdir()
                       if r.is_file() and 'rfc' in r.name and r.name.endswith('.txt'))
        if not rfc_files:
            raise FileNotFoundError('No RFC catalog files found in the data directory.')
        with resources.as_file(resources.files('vlbiplanobs.data').joinpath(sorted(rfc_files)[-1])) as rfcfile:
            return str(rfcfile)

    def _load_catalog(self) -> None:
        """Fills the instance from the (cached) parsed RFC rows, applying the flux filter.

        include_missing=True (phase cal mode): keep sources with no band measurement (negative flux),
        filter only those with a measured flux below threshold.
        include_missing=False (fringe finder mode): exclude both missing and low-flux sources.

        Raises
        ------
        FileNotFoundError
            If catalog file not found.
        RuntimeError
            If error reading the catalog file.
        """
        rows = _parse_rfc_file(self._get_catalog_path())
        unresolved = rows.flux[:, _BAND_INDEX[self._band], 1]
        if self._include_missing:
            keep = ~((unresolved >= 0) & (unresolved < self._min_flux))
        else:
            keep = ~((unresolved < 0) | (unresolved < self._min_flux))

        indices = np.flatnonzero(keep)
        # Flux row views are read-only (they view the shared cached array), so they cannot be mutated.
        sources = [CalibratorSource(rows.names[i], rows.ivsnames[i], rows.ra_deg[i], rows.dec_deg[i],
                                    int(rows.n_obs[i]), rows.flux[i, :, 0], rows.flux[i, :, 1],
                                    bool(rows.is_calibrator[i])) for i in indices]
        self._set_sources(sources, ra_arr=rows.ra_deg[indices], dec_arr=rows.dec_deg[indices])

    def _set_sources(self, sources: list[CalibratorSource], ra_arr: Optional[np.ndarray] = None,
                     dec_arr: Optional[np.ndarray] = None) -> None:
        """Sets the source list and rebuilds the name indices and coordinate arrays.

        Parameters
        ----------
        sources : list[CalibratorSource]
            Sources owned by this instance.
        ra_arr, dec_arr : np.ndarray, optional
            Pre-computed RA/Dec arrays (degrees) aligned with `sources`. Built from the sources if None.
        """
        self._sources = sources
        self._name_index = {}
        self._ivsname_index = {}
        # Keep the first occurrence for duplicated names (same behaviour as a linear search).
        for s in sources:
            self._name_index.setdefault(s.name.upper(), s)
            self._ivsname_index.setdefault(s.ivsname.upper(), s)

        if ra_arr is None or dec_arr is None:
            ra_arr = np.array([s.ra_deg for s in sources], dtype=np.float64)
            dec_arr = np.array([s.dec_deg for s in sources], dtype=np.float64)

        self._ra_arr = np.array(ra_arr, dtype=np.float64)
        self._dec_arr = np.array(dec_arr, dtype=np.float64)

    def _new_with_sources(self, sources: list[CalibratorSource]) -> Self:
        """Returns a new catalog with the same settings as this one but the given sources.

        Parameters
        ----------
        sources : list[CalibratorSource]
            Sources for the new catalog.

        Returns
        -------
        RFCCatalog
            New catalog instance (the file is not re-read).
        """
        new_catalog = object.__new__(self.__class__)
        new_catalog._min_flux = self._min_flux
        new_catalog._band = self._band
        new_catalog._include_missing = self._include_missing
        new_catalog._catalog_filename = self._catalog_filename
        new_catalog._set_sources(sources)
        return new_catalog

    @property
    def sources(self) -> list[CalibratorSource]:
        """Get the list of calibrator sources.

        Returns
        -------
        list[CalibratorSource]
            List of calibrator sources.
        """
        return self._sources

    @property
    def n_sources(self) -> int:
        """Get the number of sources in the catalog.

        Returns
        -------
        int
            Number of sources.
        """
        return len(self._sources)

    def get_source(self, name: str) -> Optional[CalibratorSource]:
        """Get a source by exact (case-insensitive) J2000 name or IVS name.

        Only sources passing this catalog's flux filter are searched. No partial/substring matching is
        done, so a non-matching name never returns a different source.

        Parameters
        ----------
        name : str
            Source name to search for.

        Returns
        -------
        CalibratorSource or None
            The matching source, or None if not found.
        """
        name_upper = name.strip().upper()
        if name_upper in self._name_index:
            return self._name_index[name_upper]

        return self._ivsname_index.get(name_upper)

    def calibrators_only(self) -> Self:
        """Return a new catalog containing only calibrator sources.

        Returns
        -------
        RFCCatalog
            New catalog with only calibrator sources.
        """
        return self._new_with_sources([s for s in self._sources if s.is_calibrator])

    def brighter_than(self, flux: float, band: Optional[str] = None) -> Self:
        """Return a new catalog with sources brighter than a threshold.

        Parameters
        ----------
        flux : float
            Minimum flux threshold in Jy.
        band : str, optional
            RFC band code to check. Default is None (use catalog default).

        Returns
        -------
        RFCCatalog
            New catalog with brighter sources.
        """
        check_band = band if band is not None else self._band
        return self._new_with_sources([s for s in self._sources if s.unresolved_flux(check_band) >= flux])

    def _get_coord_arrays(self) -> tuple[np.ndarray, np.ndarray]:
        """Get coordinate arrays for all sources.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            (ra_deg_array, dec_deg_array).
        """
        return self._ra_arr, self._dec_arr


def _angular_separation(ra1_deg: float, dec1_deg: float, ra2_arr: np.ndarray, dec2_arr: np.ndarray) -> np.ndarray:
    """Vectorized angular separation calculation using haversine formula.

    Parameters
    ----------
    ra1_deg : float
        Right ascension of first point in degrees.
    dec1_deg : float
        Declination of first point in degrees.
    ra2_arr : np.ndarray
        Array of right ascensions in degrees.
    dec2_arr : np.ndarray
        Array of declinations in degrees.

    Returns
    -------
    np.ndarray
        Angular separations in degrees.
    """
    ra1_rad, dec1_rad = np.radians(ra1_deg), np.radians(dec1_deg)
    ra2_rad, dec2_rad = np.radians(ra2_arr), np.radians(dec2_arr)

    # Vectorized haversine formula
    dlon = ra2_rad - ra1_rad
    dlat = dec2_rad - dec1_rad

    a = np.sin(dlat / 2.0) ** 2 + np.cos(dec1_rad) * np.cos(dec2_rad) * np.sin(dlon / 2.0) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(a)))


def _batch_altaz_erfa(ra_rad: np.ndarray, dec_rad: np.ndarray, times: Time,
                      station) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute elevation, azimuth, and hour angle for all sources at all times using ERFA.

    Bypasses astroplan/astropy per-source overhead by computing ERFA astrometry params
    for all time steps at once, then transforming all sources vectorized (times x sources broadcast).

    Parameters
    ----------
    ra_rad : np.ndarray
        Source right ascensions in radians.
    dec_rad : np.ndarray
        Source declinations in radians.
    times : Time
        Array of observation times.
    station : Station
        Station object with location information.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        (elevation_deg, azimuth_deg, ha_hours) each shaped (n_times, n_sources).
    """
    loc = station.location
    lon_rad, lat_rad = loc.lon.rad, loc.lat.rad
    height_m = loc.height.to(u.m).value
    utc1, utc2 = times.utc.jd1, times.utc.jd2
    dut1 = times.delta_ut1_utc
    # apco13/atciq/atioq are ufuncs: astrom has shape (n_times,), broadcast against (n_src,) -> (n_times, n_src).
    astrom, _eo = erfa.apco13(utc1, utc2, dut1, lon_rad, lat_rad, height_m, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    astrom = np.atleast_1d(astrom)[:, np.newaxis]
    ri, di = erfa.atciq(ra_rad[np.newaxis, :], dec_rad[np.newaxis, :], 0.0, 0.0, 0.0, 0.0, astrom)
    az, zen, ha, _dec, _ra = erfa.atioq(ri, di, astrom)
    elev_out = np.degrees(np.pi / 2.0 - zen)
    az_out = np.degrees(az)
    ha_out = (np.degrees(ha) / 15.0) % 24.0
    return elev_out, az_out, ha_out


def _station_observable_mask(elev: np.ndarray, az: np.ndarray, ha_hours: np.ndarray,
                             dec_deg: np.ndarray, station) -> np.ndarray:
    """Build boolean observable mask (n_times, n_sources) respecting station mount constraints.

    Uses the same constraint logic as astroplan: for ALTAZ mounts checks azimuth and elevation
    limits; for EQUAT mounts checks hour angle (in [0,24h) with wrapping), declination, and
    a minimum 5-degree elevation.

    Parameters
    ----------
    elev : np.ndarray
        Elevation angles in degrees, shape (n_times, n_sources).
    az : np.ndarray
        Azimuth angles in degrees, shape (n_times, n_sources).
    ha_hours : np.ndarray
        Hour angles in hours, shape (n_times, n_sources).
    dec_deg : np.ndarray
        Source declinations in degrees.
    station : Station
        Station object with mount constraints.

    Returns
    -------
    np.ndarray
        Boolean mask of shape (n_times, n_sources).
    """
    mount = station.mount
    if mount.mount_type == MountType.ALTAZ:
        az_min = mount.ax1.limits[0].to(u.deg).value
        az_max = mount.ax1.limits[1].to(u.deg).value
        el_min = mount.ax2.limits[0].to(u.deg).value
        el_max = mount.ax2.limits[1].to(u.deg).value
        mask = (az > az_min) & (az < az_max) & (elev > el_min) & (elev < el_max)
    else:
        ha_min = mount.ax1.limits[0].to(u.hourangle).value
        ha_max = mount.ax1.limits[1].to(u.hourangle).value
        dec_min = mount.ax2.limits[0].to(u.deg).value
        dec_max = mount.ax2.limits[1].to(u.deg).value
        ha_ok = ((ha_min < ha_hours) & (ha_hours < ha_max)) | ((ha_min + 24.0 < ha_hours) & (ha_hours < ha_max + 24.0))
        dec_ok = (dec_deg >= dec_min) & (dec_deg <= dec_max)
        mask = ha_ok & dec_ok & (elev > 5.0)

    # Apply the azimuth-dependent local horizon, if defined for this station
    if station.horizon is not None:
        mask &= elev > station.horizon_min_elevation(az)

    return mask


def get_fringe_finder_sources(
        stations: Stations, times: Time,
        min_elevation: u.Quantity = _DEFAULT_MIN_ELEVATION,
        min_flux: u.Quantity = _DEFAULT_MIN_FLUX,
        catalog: Optional[RFCCatalog] = None,
        require_all_stations: bool = True
) -> tuple[list[CalibratorSource], list[float], Optional[list[tuple[int, int, bool, float]]]]:
    """Find fringe finder sources visible by given stations during the given time range.

    Uses ERFA directly for vectorized elevation computation across all sources,
    replacing per-source astroplan calls for major speedup.

    Parameters
    ----------
    stations : Stations
        Stations object containing participating antennas.
    times : Time
        Observation times.
    min_elevation : Quantity, optional
        Minimum elevation threshold. Default is 20 degrees.
    min_flux : Quantity, optional
        Minimum flux threshold in Jy. Default is 1.0 Jy.
    catalog : RFCCatalog, optional
        Pre-loaded catalog. If None, creates a new one.
    require_all_stations : bool, optional
        If True, requires source visible by all stations. Default is True.

    Returns
    -------
    tuple[list[CalibratorSource], list[float], Optional[list[tuple[int, int, bool, float]]]]
        (sources, min_elevs, antenna_visibility) where antenna_visibility is
        None if require_all_stations=True, otherwise a list of
        (n_visible, n_total, visible_all_times, min_elev_all) tuples.
    """
    if catalog is None:
        catalog = RFCCatalog(min_flux=min_flux, band='c')

    station_list = stations.stations
    n_stations = len(station_list)
    if n_stations == 0:
        return [], [], None

    min_el_deg = float(min_elevation.to(u.deg).value) if hasattr(min_elevation, 'to') else float(min_elevation)
    valid_sources = catalog.sources
    if not valid_sources:
        return [], [], None if not require_all_stations else None

    n_sources = len(valid_sources)
    n_times = len(times)
    ra_rad = np.radians(np.array([s.ra_deg for s in valid_sources], dtype=np.float64))
    dec_rad = np.radians(np.array([s.dec_deg for s in valid_sources], dtype=np.float64))

    elev_matrices = np.empty((n_stations, n_times, n_sources))
    meets_all = np.empty((n_stations, n_times, n_sources), dtype=bool)

    dec_deg = np.degrees(dec_rad)
    for s_idx, station in enumerate(station_list):
        elev, az, ha_hours = _batch_altaz_erfa(ra_rad, dec_rad, times, station)
        obs_mask = _station_observable_mask(elev, az, ha_hours, dec_deg, station)
        elev_matrices[s_idx] = elev
        meets_all[s_idx] = obs_mask & (elev >= min_el_deg)

    if require_all_stations:
        all_times_per_station = np.all(meets_all, axis=1)
        passes_all = np.all(all_times_per_station, axis=0)
        visible_idx = np.where(passes_all)[0]

        if not visible_idx.size:
            return [], [], None

        min_elevs = np.min(elev_matrices[:, :, visible_idx], axis=(0, 1))
        sort_order = np.argsort(-min_elevs)
        sorted_idx = visible_idx[sort_order]
        return ([valid_sources[i] for i in sorted_idx],
                min_elevs[sort_order].tolist(), None)
    else:
        visible_per_station = np.any(meets_all, axis=1)
        visible_all_times_per_station = np.all(meets_all, axis=1)
        visible_counts = np.sum(visible_per_station, axis=0).astype(np.int32)
        any_visible = visible_counts > 0
        visible_idx = np.where(any_visible)[0]

        if not visible_idx.size:
            return [], [], []

        n_vis = len(visible_idx)
        result_min_elevs = np.empty(n_vis)
        result_vis_all_times = np.empty(n_vis, dtype=bool)
        result_min_elev_all = np.empty(n_vis)

        for k, src_i in enumerate(visible_idx):
            sta_mask = visible_per_station[:, src_i]
            first_sta = np.argmax(sta_mask)
            time_mask = meets_all[first_sta, :, src_i]
            result_min_elevs[k] = np.min(elev_matrices[first_sta, time_mask, src_i])
            result_vis_all_times[k] = np.all(visible_all_times_per_station[sta_mask, src_i])
            result_min_elev_all[k] = np.min(elev_matrices[sta_mask, :, src_i])

        vc = visible_counts[visible_idx]
        categories = np.where((vc == n_stations) & result_vis_all_times, 0,
                     np.where(vc == n_stations, 1, 2))
        sort_order = np.lexsort((-result_min_elev_all, -vc, categories))

        sorted_idx = visible_idx[sort_order]
        sorted_sources = [valid_sources[i] for i in sorted_idx]
        sorted_min_elevs = result_min_elevs[sort_order].tolist()
        sorted_antenna_vis = [(int(vc[sort_order[j]]), n_stations,
                               bool(result_vis_all_times[sort_order[j]]),
                               float(result_min_elev_all[sort_order[j]])) for j in range(n_vis)]

        return sorted_sources, sorted_min_elevs, sorted_antenna_vis


def get_nearby_sources(source: CalibratorSource | Source, max_separation: u.Quantity = 5.0 * u.deg,
                       catalog: Optional[RFCCatalog] = None,
                       n_sources: Optional[int] = None) -> list[tuple[CalibratorSource, float]]:
    """Find calibrator sources near a target source.

    Parameters
    ----------
    source : CalibratorSource or Source
        Target source.
    max_separation : Quantity, optional
        Maximum angular separation. Default is 5 degrees.
    catalog : RFCCatalog, optional
        Pre-loaded catalog. If None, creates a new one.
    n_sources : int, optional
        Maximum number of sources to return. If None, returns all.

    Returns
    -------
    list[tuple[CalibratorSource, float]]
        List of (source, separation_deg) tuples sorted by separation.
    """
    if catalog is None:
        catalog = RFCCatalog()

    max_sep_deg = max_separation.to(u.deg).value if hasattr(max_separation, 'to') else max_separation
    ra_arr, dec_arr = catalog._get_coord_arrays()
    if not ra_arr.size:
        return []

    separations_deg = _angular_separation(source.coord.ra.deg, source.coord.dec.deg, ra_arr, dec_arr)
    valid_mask = (separations_deg <= max_sep_deg) & (separations_deg > 0)

    valid_indices = np.where(valid_mask)[0]
    nearby = [(catalog.sources[i], separations_deg[i]) for i in valid_indices]
    nearby.sort(key=lambda x: x[1])
    return nearby[:n_sources] if n_sources is not None else nearby


def select_phase_calibrator(target: Source, band: str, max_separation: u.Quantity = 5.0 * u.deg,
                            catalog: Optional[RFCCatalog] = None) -> Optional[CalibratorSource]:
    """Automatically select the best phase calibrator for a target source.

    Prioritises the brightest unresolved emission among the closest sources.
    The scoring formula is ``score = unresolved_flux / (1 + separation_deg)``
    so that a bright, nearby, compact source is preferred.

    Parameters
    ----------
    target : Source
        The target source to find a phase calibrator for.
    band : str
        Observing band in RFC letter code (s/c/x/u/k) or wavelength string (e.g. '6cm').
    max_separation : Quantity
        Maximum angular separation from target (default 5 deg).
    catalog : RFCCatalog or None
        Pre-loaded catalog.  If None a new one is created with min_flux=0.

    Returns
    -------
    CalibratorSource or None
        The best phase calibrator, or None if no candidates exist.
    """
    rfc_band = _wavelength_to_rfc_band(band)
    if catalog is None:
        catalog = RFCCatalog(min_flux=0.0 * u.Jy, band=rfc_band)

    nearby = get_nearby_sources(target, max_separation=max_separation, catalog=catalog)
    if not nearby:
        return None

    target_names = {target.name.upper()}
    if hasattr(target, 'other_names') and target.other_names:
        target_names.update(n.upper() for n in target.other_names)

    best_src, best_score = None, -1.0
    for src, sep_deg in nearby:
        # Skip the target source itself
        if src.name.upper() in target_names or src.ivsname.upper() in target_names:
            continue
        flux_unres = src.unresolved_flux(rfc_band)
        if flux_unres <= 0:
            _, flux_unres = src.get_flux_at_band(band)
        if flux_unres <= 0:
            continue
        # Compactness bonus: ratio of unresolved to resolved flux (capped at 1)
        flux_res = src.resolved_flux(rfc_band)
        compactness = min(flux_unres / flux_res, 1.0) if flux_res > 0 else 1.0
        score = flux_unres * compactness / (1.0 + sep_deg)
        if score > best_score:
            best_score, best_src = score, src

    return best_src


def select_check_source(target: Source, phase_cal: Source, band: str,
                        max_separation: u.Quantity = 5.0 * u.deg,
                        catalog: Optional[RFCCatalog] = None) -> Optional[CalibratorSource]:
    """Automatically select the best check source for a target/phase-cal pair.

    Prefers a source that is:
    - close to the target,
    - at roughly the same distance from the target as the phase calibrator,
    - on the same side of the sky (similar position angle),
    - compact (high unresolved / resolved ratio).
    It can be weaker than the phase calibrator.

    Parameters
    ----------
    target : Source
        The target source.
    phase_cal : Source
        The already-selected phase calibrator.
    band : str
        Observing band (RFC letter or wavelength string).
    max_separation : Quantity
        Maximum angular separation from target.
    catalog : RFCCatalog or None
        Pre-loaded catalog.

    Returns
    -------
    CalibratorSource or None
        The best check source, or None if no candidates exist.
    """
    rfc_band = _wavelength_to_rfc_band(band)
    if catalog is None:
        catalog = RFCCatalog(min_flux=0.0 * u.Jy, band=rfc_band)

    nearby = get_nearby_sources(target, max_separation=max_separation, catalog=catalog)
    if not nearby:
        return None

    # Reference geometry: target → phase_cal
    pc_sep = float(target.coord.separation(phase_cal.coord).deg)
    pc_pa = float(target.coord.position_angle(phase_cal.coord).deg)

    # Build set of names to exclude (target + phase cal)
    exclude_names: set[str] = {target.name.upper(), phase_cal.name.upper()}
    if hasattr(target, 'other_names') and target.other_names:
        exclude_names.update(n.upper() for n in target.other_names)
    if hasattr(phase_cal, 'other_names') and phase_cal.other_names:
        exclude_names.update(n.upper() for n in phase_cal.other_names)

    best_src, best_score = None, -1.0
    for src, sep_deg in nearby:
        # Exclude target and phase calibrator
        if src.name.upper() in exclude_names:
            continue
        if hasattr(src, 'ivsname') and src.ivsname.upper() in exclude_names:
            continue
        flux_unres = src.unresolved_flux(rfc_band)
        if flux_unres <= 0:
            _, flux_unres = src.get_flux_at_band(band)
        if flux_unres <= 0:
            continue

        flux_res = src.resolved_flux(rfc_band)
        compactness = min(flux_unres / flux_res, 1.0) if flux_res > 0 else 1.0

        # Distance-match bonus: prefer sources at similar distance as phase cal
        dist_match = 1.0 / (1.0 + abs(sep_deg - pc_sep))

        # Direction bonus: prefer sources on the same side as the phase cal
        src_pa = float(target.coord.position_angle(src.coord).deg)
        pa_diff = abs(src_pa - pc_pa) % 360
        if pa_diff > 180:
            pa_diff = 360 - pa_diff
        direction_bonus = 1.0 / (1.0 + pa_diff / 90.0)

        # Combined score: compactness matters most, proximity next, geometry bonus
        score = compactness * (0.3 * flux_unres + 0.7) * dist_match * direction_bonus / (1.0 + sep_deg)
        if score > best_score:
            best_score, best_src = score, src

    return best_src


def _wavelength_to_rfc_band(band: str) -> str:
    """Convert a wavelength string like '6cm' to an RFC band letter like 'c'.

    Also accepts bare RFC band letters ('s', 'c', 'x', 'u', 'k') unchanged.

    Parameters
    ----------
    band : str
        Wavelength string (e.g. '6cm', '18cm') or RFC letter.

    Returns
    -------
    str
        Single-letter RFC band code.
    """
    if band in _BAND_INDEX:
        return band
    if band in _RFC_BANDS:
        return _RFC_BANDS[band]
    if band in _WAVELENGTH_BANDS:
        return _WAVELENGTH_BANDS[band]
    if band.lower().endswith('cm'):
        rounded = _round_to_nearest_wavelength(band)
        if rounded in _WAVELENGTH_BANDS:
            return _WAVELENGTH_BANDS[rounded]
    return 'c'


# Single source of truth for the CLI defaults of the fringe-finder and phase-calibrator searches.
# Used by the 'planobs fringefinders/phasecals' option handlers in cli.py (also behind main_fringe/main_phasecal).
FRINGE_DEFAULT_MIN_FLUX_JY: float = 0.5
FRINGE_DEFAULT_MIN_ELEVATION_DEG: float = 20.0
FRINGE_DEFAULT_MAX_LINES: int = 20
PHASECAL_DEFAULT_MAX_SEPARATION_DEG: float = 5.0
PHASECAL_DEFAULT_MIN_FLUX_JY: float = 0.1


def _print_error(message: str, as_json: bool) -> None:
    """Prints an error message either as a JSON object ({"error": message}) or as bold red rich text.

    Inputs
        message : str — the error message.
        as_json : bool — if True, print JSON to stdout; otherwise print rich-formatted text.
    """
    if as_json:
        print(json.dumps({"error": message}, indent=2))
    else:
        rprint(f"[bold red]{message}[/bold red]")


def _display_fluxes(src: Source, band: Optional[str]) -> tuple[float, float]:
    """Returns the (total, unresolved) flux in Jy to display for a source.

    Inputs
        src : Source — a calibrator source with flux information.
        band : Optional[str] — observing band; if None, the maximum flux over all bands is used.

    Returns
        tuple[float, float] — (total_flux_jy, unresolved_flux_jy); 0.0 when no positive flux is available.
    """
    if band:
        return src.get_flux_at_band(band)

    total_flux = float(np.max(src.flux_resolved)) if np.any(src.flux_resolved > 0) else 0.0
    unresolved_flux = float(np.max(src.flux_unresolved)) if np.any(src.flux_unresolved > 0) else 0.0
    return total_flux, unresolved_flux


def _find_station_codename(name: str) -> Optional[str]:
    """Resolves an antenna name or codename to its codename in the loaded station catalog (obs._STATIONS).

    Inputs
        name : str — antenna codename or full name. Exact match first, then case-insensitive
                     match against codenames and then against full names.

    Returns
        Optional[str] — the station codename, or None if the antenna is unknown.
    """
    try:
        return obs._STATIONS[name.strip()].codename
    except KeyError:
        pass

    name_upper = name.strip().upper()
    for codename in obs._STATIONS.station_codenames:
        if codename.upper() == name_upper:
            return obs._STATIONS[codename].codename

    for full_name in obs._STATIONS.station_names:
        if full_name.upper() == name_upper:
            return obs._STATIONS[full_name].codename

    return None


def _visibility_text(visibility: tuple, long_form: bool) -> str:
    """Formats the per-source antenna visibility summary.

    Inputs
        visibility : tuple — (visible_count, total_count, visible_all_times, min_elev_all) as returned
                     by get_fringe_finder_sources.
        long_form : bool — True for the table wording ("all, all time"), False for the JSON wording
                    ("all ant. all time").

    Returns
        str — human-readable visibility summary.
    """
    visible_count, total_count, visible_all_times, _ = visibility
    if visible_count == total_count and visible_all_times:
        return "all, all time" if long_form else "all ant. all time"
    if visible_count == total_count:
        return "all, partial time" if long_form else "all ant. partial time"
    return f"{visible_count}/{total_count} antennas" if long_form else f"{visible_count}/{total_count} ant."


def run_fringe_finders(*, starttime: str | Time, duration: float, networks: Optional[list[str]] = None,
                       stations: Optional[list[str]] = None, min_flux: float = FRINGE_DEFAULT_MIN_FLUX_JY,
                       min_elevation: float = FRINGE_DEFAULT_MIN_ELEVATION_DEG,
                       max_lines: int = FRINGE_DEFAULT_MAX_LINES, require_all: bool = False,
                       band: Optional[str] = None, station_catalog: Optional[str] = None,
                       as_json: bool = False) -> int:
    """Searches for fringe-finder candidates for the given antennas/time range and prints the results.

    Side effect: replaces the global station catalog (obs._STATIONS) with the one from `station_catalog`
    (or the default catalog if None).

    Inputs
        starttime : str | Time — start of the observation (UTC), e.g. '2025-03-15 08:00'.
        duration : float — duration of the observation in hours.
        networks : Optional[list[str]] — VLBI network names; their default stations are included.
        stations : Optional[list[str]] — antenna codenames or names to include.
        min_flux : float — minimum unresolved flux in Jy.
        min_elevation : float — minimum elevation in degrees.
        max_lines : int — maximum number of sources to print.
        require_all : bool — require the source to be visible by all antennas.
        band : Optional[str] — band for the flux display (e.g. '6cm'); None shows the maximum over all bands.
        station_catalog : Optional[str] — path to a custom station catalog file.
        as_json : bool — print JSON instead of a rich table.

    Returns
        int — process exit code: 0 on success (also when no candidates are found), 1 on input errors.
    """
    if networks is None and stations is None:
        _print_error("You need to provide at least a VLBI network or a list of antennas.", as_json)
        return 1

    obs._STATIONS = obs.Stations(filename=station_catalog)
    stations_list: list[str] = []
    unknown_networks = [n for n in networks or [] if n not in obs._NETWORKS]
    if unknown_networks:
        n_networks = len(unknown_networks)
        _print_error(f"The network{'s' if n_networks > 1 else ''} {', '.join(unknown_networks)} "
                     f"{'are' if n_networks > 1 else 'is'} not known.", as_json)
        return 1

    for n in networks or []:
        for s in obs._NETWORKS[n].station_codenames:
            if s not in stations_list:
                stations_list.append(s)

    for s in stations or []:
        a_station = _find_station_codename(s)
        if a_station is None:
            _print_error(f"The station {s} is not known.", as_json)
            return 1
        if a_station not in stations_list:
            stations_list.append(a_station)

    stations_obj = obs._STATIONS.filter_antennas(stations_list)
    if not stations_obj:
        _print_error("No valid antennas have been selected.", as_json)
        return 1

    _log.info("Fringe-finder search: stations=%s start=%s duration=%sh min_flux=%sJy min_elev=%sdeg",
              stations_list, starttime, duration, min_flux, min_elevation)
    times = Time(starttime, scale='utc') + np.arange(0, duration + 0.1, 0.1) * u.hour
    sources, min_elevs, antenna_visibility = get_fringe_finder_sources(stations_obj, times,
                                                                       min_elevation=min_elevation * u.deg,
                                                                       min_flux=min_flux * u.Jy,
                                                                       require_all_stations=require_all)
    if not sources:
        _print_error(f"No fringe finder candidates found above {min_elevation} degrees elevation "
                     f"and with a unresolved flux above {min_flux} Jy.", as_json)
        return 0

    max_display = min(max_lines, len(sources))
    if as_json:
        result_data = []
        for i in range(max_display):
            src, min_elev = sources[i], min_elevs[i]
            total_flux, unresolved_flux = _display_fluxes(src, band)
            entry = {"name": src.name, "ivs_name": src.ivsname,
                     "min_elevation_deg": min_elev if min_elev > 0.0 else 0,
                     "total_flux_jy": total_flux, "unresolved_flux_jy": unresolved_flux,
                     "bands": src.get_observed_bands(), "astrogeo_url": src.get_astrogeo_link()}
            if antenna_visibility is not None:
                entry["antenna_visibility"] = _visibility_text(antenna_visibility[i], long_form=False)
            result_data.append(entry)

        result = {"min_elevation_deg": min_elevation, "min_flux_jy": min_flux,
                  "require_all_stations": require_all, "sources": result_data,
                  "total_found": len(sources), "shown": max_display}
        print(json.dumps(result, indent=2))
        return 0

    rprint(f"\n[bold green]Found {len(sources)} fringe finder candidates above "
           f"{min_elevation} degrees elevation and with a unresolved flux "
           f"above {min_flux} Jy:[/bold green]")
    table = Table(show_header=True, header_style="bold", show_lines=False, box=box.SIMPLE)
    table.add_column("Name", style="", width=17)
    table.add_column("IVS Name", style="", width=10)
    table.add_column("Min elev. (deg)", justify="right", style="", width=10)
    table.add_column("Total flux (Jy)", justify="right", style="", width=12)
    table.add_column("Unresolved (Jy)", justify="right", style="", width=13)
    table.add_column("Bands", justify="right", style="", width=10)
    table.add_column("url", style="", width=10)
    if antenna_visibility is not None:
        table.add_column("Antenna Visibility", justify="center", style="", width=15)

    for i in range(max_display):
        src, min_elev = sources[i], min_elevs[i]
        total_flux, unresolved_flux = _display_fluxes(src, band)
        row = [src.name, src.ivsname, f"{min_elev if min_elev > 0.0 else 0:>6.1f}",
               f"{total_flux:>8.2f}" if total_flux > 0 else "N/A",
               f"{unresolved_flux:>8.2f}" if unresolved_flux > 0 else "N/A",
               src.get_observed_bands(), f"[link={src.get_astrogeo_link()}]AstroGeo[/link]"]
        if antenna_visibility is not None:
            row.append(_visibility_text(antenna_visibility[i], long_form=True))
        table.add_row(*row)

    rprint(table)
    if len(sources) > max_display:
        rprint(f"\n... and {len(sources) - max_display} more sources.")
    return 0


def main_fringe():
    """Console entry point 'planobs_fringefinder': same options and behaviour as 'planobs fringefinders'.

    Builds the parser with cli.add_fringe_finder_arguments (single definition of the option names) and
    exits through cli.handle_fringe_finder_command with run_fringe_finders' exit code.
    """
    from vlbiplanobs import cli
    parser = argparse.ArgumentParser(description="Find fringe finder sources for VLBI observations",
                                     prog="planobs_fringefinder", formatter_class=RawTextRichHelpFormatter)
    cli.add_fringe_finder_arguments(parser)
    cli.handle_fringe_finder_command(parser.parse_args())


def _target_from_personal_catalog(catalog_file: str, name: str) -> Optional[Source]:
    """Looks up `name` in a personal source catalog (toml file, as used by 'planobs observe -sc').

    Inputs
        catalog_file : str — path to the personal source catalog (toml).
        name : str — block name or source name to look up.

    Returns
        Optional[Source] — the (first) target source of the matching block, the source with
        the given name, or None if not found.

    Raises
        FileNotFoundError — if catalog_file does not exist.
    """
    catalog = SourceCatalog(catalog_file)
    if name in catalog.blocknames:
        block = catalog[name]
        block_targets = block.sources(SourceType.TARGET) or block.sources()
        return block_targets[0] if block_targets else None

    return catalog.sources(include_calibrators=True).get(name)


def run_phasecals(*, target: str, max_separation: float = PHASECAL_DEFAULT_MAX_SEPARATION_DEG,
                  min_flux: float = PHASECAL_DEFAULT_MIN_FLUX_JY, n_sources: Optional[int] = None,
                  band: Optional[str] = None, catalog_file: Optional[str] = None,
                  source_catalog: Optional[str] = None, as_json: bool = False) -> int:
    """Searches for phase-calibrator candidates near a target source and prints the results.

    The target is resolved in this order: personal `source_catalog` (block or source name), the RFC
    catalog, and finally Source.source_from_str (coordinates or online name lookup).

    Inputs
        target : str — target source name (J2000/IVS name, block/source name in `source_catalog`, or coordinates).
        max_separation : float — maximum angular separation in degrees.
        min_flux : float — minimum unresolved flux in Jy.
        n_sources : Optional[int] — maximum number of sources to return; None returns all.
        band : Optional[str] — band for the flux display (e.g. '6cm'); None shows the maximum over all bands.
        catalog_file : Optional[str] — path to a custom RFC catalog file.
        source_catalog : Optional[str] — path to a personal source catalog (toml).
        as_json : bool — print JSON instead of a rich table.

    Returns
        int — process exit code: 0 on success, 1 if the target cannot be resolved or no candidates are found.
    """
    catalog = RFCCatalog(catalog_filename=catalog_file, band='c', min_flux=min_flux, include_missing=True)
    target_src = None
    if source_catalog is not None:
        try:
            target_src = _target_from_personal_catalog(source_catalog, target)
        except FileNotFoundError:
            _print_error(f"Source catalog file not found: {source_catalog}", as_json)
            return 1
        if target_src is None:
            rprint(f"[yellow]'{target}' not found in {source_catalog} — "
                   "falling back to the RFC catalog/online lookup.[/yellow]")

    if not target_src:
        target_src = catalog.get_source(target)
    if not target_src:
        try:
            target_src = Source.source_from_str(target, source_type=SourceType.TARGET)
        except Exception:
            _print_error(f"Target source '{target}' could not be parsed or found in the catalogs.", as_json)
            return 1

    _log.info("Phase-calibrator search: target=%s max_sep=%sdeg min_flux=%sJy n_sources=%s",
              target_src.name, max_separation, min_flux, n_sources)
    nearby = get_nearby_sources(target_src, max_separation=max_separation * u.deg, catalog=catalog,
                                n_sources=n_sources)
    target_coords = target_src.coord.to_string('hmsdms')
    if not nearby:
        _print_error(f"No phase calibrator candidates found near {target_src.name} ({target_coords}).", as_json)
        return 1

    if as_json:
        result_data = []
        for src, sep in nearby:
            total_flux, unresolved_flux = _display_fluxes(src, band)
            result_data.append({"name": src.name, "ivs_name": src.ivsname, "separation_deg": sep,
                                "total_flux_jy": total_flux, "unresolved_flux_jy": unresolved_flux,
                                "bands": src.get_observed_bands(), "astrogeo_url": src.get_astrogeo_link()})

        result = {"target_name": target_src.name, "target_coordinates": target_coords,
                  "max_separation_deg": max_separation, "min_flux_jy": min_flux,
                  "sources": result_data, "total_found": len(nearby)}
        print(json.dumps(result, indent=2))
        return 0

    rprint(f"\n[bold green]Found {len(nearby)} phase calibrator candidates near {target_src.name} "
           f"({target_coords}):[/bold green]")
    table = Table(show_header=True, header_style="bold", show_lines=False, box=box.SIMPLE)
    table.add_column("Name", style="", width=17)
    table.add_column("IVS Name", style="", width=10)
    table.add_column("Separation (deg)", justify="right", style="", width=12)
    table.add_column("Total flux (Jy)", justify="right", style="", width=12)
    table.add_column("Unresolved (Jy)", justify="right", style="", width=13)
    table.add_column("Bands", justify="right", style="", width=10)
    table.add_column("url", style="", width=10)
    for src, sep in nearby:
        total_flux, unresolved_flux = _display_fluxes(src, band)
        table.add_row(src.name, src.ivsname, f"{sep:.2f}",
                      f"{total_flux:>8.2f}" if total_flux > 0 else "N/A",
                      f"{unresolved_flux:>8.2f}" if unresolved_flux > 0 else "N/A",
                      src.get_observed_bands(), f"[link={src.get_astrogeo_link()}]AstroGeo[/link]")

    rprint(table)
    return 0


def main_phasecal():
    """Console entry point 'planobs_phasecal': same options and behaviour as 'planobs phasecals'.

    Builds the parser with cli.add_phase_cal_arguments (single definition of the option names) and
    exits through cli.handle_phase_cal_command with run_phasecals' exit code.
    """
    from vlbiplanobs import cli
    parser = argparse.ArgumentParser(description="Find phase calibrator sources near a target source",
                                     prog="planobs_phasecal", formatter_class=RawTextRichHelpFormatter)
    cli.add_phase_cal_arguments(parser)
    cli.handle_phase_cal_command(parser.parse_args())
