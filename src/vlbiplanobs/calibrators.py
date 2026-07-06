import sys
import argparse
import json
from typing import Optional, Self
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

"""Calibrator source module for finding fringe finders (bright, compact sources used to check
the fringe fit) and nearby phase calibrators for a target source, based on the RFC (Radio
Fundamental Catalog).
"""

_RFC_BANDS: dict[str, str] = {'l': 's', 's': 's', 'c': 'c', 'm': 'c', 'x': 'x', 'u': 'u', 'k': 'k', 'q': 'k'}
_DEFAULT_MIN_ELEVATION: u.Quantity = 20 * u.deg
_DEFAULT_MIN_FLUX: u.Quantity = 1.0 * u.Jy
_BAND_INDEX: dict[str, int] = {'s': 0, 'c': 1, 'x': 2, 'u': 3, 'k': 4}
_WAVELENGTH_BANDS: dict[str, str] = {'18cm': 's', '21cm': 's', '13cm': 'c', '6cm': 'c', '5cm': 'c',
                       '3.6cm': 'x', '2cm': 'u', '1.3cm': 'k', '0.7cm': 'k'}


def _round_to_nearest_wavelength(band_str: str) -> str:
    """Round an unknown wavelength string to the nearest known wavelength key in _WAVELENGTH_BANDS.

    Inputs
    - band_str : str
        Wavelength string (e.g. '6.2cm'). If it already matches a known key, or it does not
        end in 'cm', it is returned unchanged.

    Returns
    - str
        The closest known wavelength key (e.g. '6cm'), or `band_str` unchanged if it cannot
        be matched.
    """
    if band_str in _WAVELENGTH_BANDS:
        return band_str
    if not band_str.lower().endswith('cm'):
        return band_str
    return min(_WAVELENGTH_BANDS.keys(), key=lambda x: abs(float(x[:-2]) - float(band_str[:-2])))


class CalibratorSource(Source):
    """Represents a calibrator source from the RFC catalog.
    Stores the RFC-specific metadata (IVS name, number of VLBI observations, resolved and
    unresolved flux density per band, and whether it is flagged as a calibrator) on top of the
    base `Source` attributes (name, coordinates).
    """
    __slots__ = ('ivsname', 'n_observations', 'flux_resolved', 'flux_unresolved', 'is_calibrator')

    def __init__(self, name: str, ivsname: str, ra_deg: float, dec_deg: float, n_observations: int,
                 flux_resolved: np.ndarray, flux_unresolved: np.ndarray, is_calibrator: bool):
        """Initializes a calibrator source.

        Inputs
        - name : str
            Primary (IERS/RFC) name of the source.
        - ivsname : str
            IVS name of the source, as used in VLBI scheduling.
        - ra_deg : float
            Right ascension in degrees (J2000).
        - dec_deg : float
            Declination in degrees (J2000).
        - n_observations : int
            Number of VLBI observations of this source in the RFC catalog.
        - flux_resolved : numpy.ndarray
            Resolved flux density (Jy) at each of the five RFC bands (s, c, x, u, k), following
            the order in `_BAND_INDEX`.
        - flux_unresolved : numpy.ndarray
            Unresolved (core) flux density (Jy) at each of the five RFC bands, same order as
            `flux_resolved`.
        - is_calibrator : bool
            True if the source is flagged in the RFC catalog as a good calibrator ('C' entries).
        """
        super().__init__(name=name, coordinates=coord.SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg),
                     source_type=SourceType.PHASECAL)
        self.ivsname: str = ivsname
        self.n_observations: int = n_observations
        self.flux_resolved: np.ndarray = flux_resolved
        self.flux_unresolved: np.ndarray = flux_unresolved
        self.is_calibrator: bool = is_calibrator

    @property
    def ra_deg(self) -> float:
        """Right ascension of the source in degrees (J2000)."""
        return float(self.coord.ra.deg)

    @property
    def dec_deg(self) -> float:
        """Declination of the source in degrees (J2000)."""
        return float(self.coord.dec.deg)

    def unresolved_flux(self, band: str) -> float:
        """Returns the unresolved (core) flux density of the source at the given RFC band.

        Inputs
        - band : str
            RFC band letter code ('s', 'c', 'x', 'u', or 'k').

        Returns
        - float
            Unresolved flux density in Jy, or 0.0 if the band is not recognized.
        """
        idx = _BAND_INDEX.get(band)
        return float(self.flux_unresolved[idx]) if idx is not None else 0.0

    def resolved_flux(self, band: str) -> float:
        """Returns the resolved flux density of the source at the given RFC band.

        Inputs
        - band : str
            RFC band letter code ('s', 'c', 'x', 'u', or 'k').

        Returns
        - float
            Resolved flux density in Jy, or 0.0 if the band is not recognized.
        """
        idx = _BAND_INDEX.get(band)
        return float(self.flux_resolved[idx]) if idx is not None else 0.0

    def get_flux_at_band(self, band: str) -> tuple[float, float]:
        """Get resolved and unresolved flux at a specific band, interpolating if not measured.

        Accepts either a wavelength string (e.g. '6cm') or an RFC band letter. If the band has
        no direct measurement (flux is zero or missing), falls back to `_interpolate_flux`.

        Inputs
        - band : str
            Wavelength string (e.g. '18cm', '6cm') or RFC band letter ('s', 'c', 'x', 'u', 'k').

        Returns
        - tuple[float, float]
            (resolved_flux_jy, unresolved_flux_jy). Both 0.0 if the band cannot be matched or
            no flux data is available at all.
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
        """Interpolate resolved/unresolved flux for a band lacking a direct measurement.

        Uses log-log linear interpolation in wavelength between the two closest bands (by
        wavelength) that do have a positive measured unresolved flux. Falls back to a plain
        average of the two closest bands when interpolation is not possible (identical
        wavelengths, or a zero/negative flux at either band).

        Inputs
        - target_band : str
            RFC band letter code ('s', 'c', 'x', 'u', or 'k') to interpolate flux for.

        Returns
        - tuple[float, float]
            (resolved_flux_jy, unresolved_flux_jy) interpolated at the target band. Both 0.0
            if no other band has measured flux data, or the target band's representative
            wavelength cannot be determined.
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

    def get_skycoord(self) -> coord.SkyCoord:
        """Returns the source coordinates as an astropy.coordinates.SkyCoord object."""
        return self.coord

    def get_astrogeo_link(self) -> str:
        """Builds the AstroGeo calibrator-search URL for this source's coordinates.

        Returns
        - str
            URL to the AstroGeo (astrogeo.org) calibrator search form, pre-filled with this
            source's RA/Dec.
        """
        # Handle negative declinations properly - only degrees should be negative
        ra_h, ra_m, ra_s = self.coord.ra.hms
        dec_d, dec_m, dec_s = self.coord.dec.dms

        # For negative declinations, make minutes and seconds positive
        if dec_d < 0:
            dec_m = abs(dec_m)
            dec_s = abs(dec_s)

        fmt = "ra={:02.0f}:{:02.0f}:{:06.3f}&dec={:+03.0f}:{:02.0f}:{:06.3f}&num_sou=1&format=html"
        source_coord_str = parse.quote(fmt.format(ra_h, ra_m, ra_s, dec_d, dec_m, dec_s), safe='=&')
        return f"http://astrogeo.org/cgi-bin/calib_search_form.csh?{source_coord_str}"

    def get_observed_bands(self) -> str:
        """Returns a comma-separated, upper-case list of RFC bands with measured flux data.

        Returns
        - str
            Comma-separated band letters (e.g. 'S,C,X') for which either resolved or
            unresolved flux is positive. Empty string if no band has data.
        """
        bands = []
        for band in ['s', 'c', 'x', 'u', 'k']:
            idx = _BAND_INDEX[band]
            if self.flux_unresolved[idx] > 0 or self.flux_resolved[idx] > 0:
                bands.append(band.upper())
        return ','.join(bands)


class RFCCatalog:
    """RFC (Radio Fundamental Catalog) of VLBI calibrator sources.
    Loads and parses the catalog text file into a list of `CalibratorSource` objects, filtered
    by a minimum unresolved flux threshold at a given band, and provides lookup and filtering
    helpers (by name, by flux, calibrators-only).
    """
    __slots__ = ('_sources', '_min_flux', '_band', '_catalog_filename',
                 '_name_index', '_ivsname_index', '_ra_arr', '_dec_arr', '_include_missing')

    def __init__(self, catalog_filename: Optional[str] = None, min_flux: u.Quantity = _DEFAULT_MIN_FLUX,
                 band: str = 'c', include_missing: bool = False):
        """Initializes the catalog by loading and parsing the RFC catalog file.

        Inputs
        - catalog_filename : str [OPTIONAL]
            Path to a custom RFC catalog text file. If not provided, the most recent
            'rfc*.txt' file bundled in the `vlbiplanobs.data` package is used.
        - min_flux : astropy.units.Quantity [OPTIONAL]
            Minimum unresolved flux density (default in Jy) a source must have at `band` to be
            kept in the catalog. Defaults to 1.0 Jy.
        - band : str [OPTIONAL]
            RFC band letter code ('s', 'c', 'x', 'u', or 'k') used to apply the `min_flux` cut.
            Defaults to 'c'.
        - include_missing : bool [OPTIONAL]
            If True, sources with no flux measurement at `band` are kept (only sources with a
            measured flux below `min_flux` are excluded); used for phase-calibrator mode. If
            False (default), sources with no measurement or with flux below `min_flux` are both
            excluded; used for fringe-finder mode.

        Raises
        - FileNotFoundError: If no RFC catalog file can be found, or the given/default file
          does not exist.
        - RuntimeError: If the catalog file exists but cannot be parsed.
        """
        self._sources: list[CalibratorSource] = []
        self._min_flux = min_flux.to(u.Jy).value if hasattr(min_flux, 'to') else min_flux
        self._band = band
        self._include_missing = include_missing
        self._catalog_filename = catalog_filename
        self._name_index: dict[str, CalibratorSource] = {}
        self._ivsname_index: dict[str, CalibratorSource] = {}
        self._ra_arr: np.ndarray = np.array([], dtype=np.float64)
        self._dec_arr: np.ndarray = np.array([], dtype=np.float64)
        self._load_catalog()

    def _get_catalog_path(self) -> str:
        """Resolves the path to the RFC catalog file to load.

        Returns
        - str
            `self._catalog_filename` if set; otherwise the path to the most recent
            'rfc*.txt' file bundled in the `vlbiplanobs.data` package (sorted by filename).

        Raises
        - FileNotFoundError: If no bundled RFC catalog file is found.
        """
        if self._catalog_filename is not None:
            return self._catalog_filename
        rfc_files = tuple(r.name for r in resources.files('vlbiplanobs.data').iterdir()
                       if r.is_file() and 'rfc' in r.name and r.name.endswith('.txt'))
        if len(rfc_files) == 0:
            raise FileNotFoundError('No RFC catalog files found in the data directory.')
        with resources.as_file(resources.files('vlbiplanobs.data').joinpath(sorted(rfc_files)[-1])) as rfcfile:
            return str(rfcfile)

    def _load_catalog(self) -> None:
        """Reads the RFC catalog file, parses it, and populates `self._sources` and the
        name/IVS-name/coordinate lookup indices, applying the `_min_flux`/`_band`/
        `_include_missing` filters from `__init__`.

        Raises
        - FileNotFoundError: If the catalog file does not exist.
        - RuntimeError: If a line-parsing error occurs while reading the file.
        """
        catalog_path = self._get_catalog_path()
        try:
            band_idx = _BAND_INDEX[self._band]
            with open(catalog_path, 'rt') as fin:
                lines = fin.readlines()

            valid_lines = [line for line in lines
                          if line and line[0] not in '#U' and len(line) >= 100]

            parsed_data = []
            for line in valid_lines:
                cols = line.split()
                if len(cols) < 25:
                    continue

                try:
                    flux_values = []
                    for f_res, f_unres in zip(cols[13:23:2], cols[14:24:2]):
                        flux_res = 0.0 if f_res[0] == '<' else float(f_res)
                        flux_unres = 0.0 if f_unres[0] == '<' else float(f_unres)
                        flux_values.extend([flux_res, flux_unres])

                    flux_array = np.array(flux_values, dtype=np.float32).reshape(5, 2)
                    # Skip sources with flux below threshold.
                    # include_missing=True (phase cal mode): keep sources with no band measurement,
                    # filter only those with a measured flux below threshold.
                    # include_missing=False (fringe finder mode): exclude both missing and low-flux sources.
                    if self._include_missing:
                        if flux_array[band_idx, 1] >= 0 and flux_array[band_idx, 1] < self._min_flux:
                            continue
                    else:
                        if flux_array[band_idx, 1] < 0 or flux_array[band_idx, 1] < self._min_flux:
                            continue
                except (ValueError, IndexError):
                    continue

                try:
                    ra_deg = (float(cols[3]) + float(cols[4]) / 60.0 + float(cols[5]) / 3600.0) * 15.0
                    dec_deg = (1.0 if cols[6][0] != '-' else -1.0) * \
                        (abs(float(cols[6])) + float(cols[7]) / 60.0 + float(cols[8]) / 3600.0)
                except (ValueError, IndexError):
                    continue

                n_obs = int(cols[12]) if cols[12].lstrip('-').isdigit() else 0
                parsed_data.append((cols[2], cols[1], ra_deg, dec_deg, n_obs,
                                    flux_array[:, 0], flux_array[:, 1], cols[0] == 'C'))

            self._sources = [CalibratorSource(name, ivs, ra, dec, nobs, flux_r, flux_u, is_cal)
                             for name, ivs, ra, dec, nobs, flux_r, flux_u, is_cal in parsed_data]

            self._name_index = {s.name: s for s in self._sources}
            self._ivsname_index = {s.ivsname: s for s in self._sources}

            if self._sources:
                coords = np.array([(s.ra_deg, s.dec_deg) for s in self._sources], dtype=np.float64)
                self._ra_arr = coords[:, 0].copy()
                self._dec_arr = coords[:, 1].copy()
            else:
                self._ra_arr = np.array([], dtype=np.float64)
                self._dec_arr = np.array([], dtype=np.float64)
        except FileNotFoundError:
            raise FileNotFoundError(f'RFC catalog file not found: {catalog_path}')
        except Exception as e:
            raise RuntimeError(f'Error parsing RFC catalog file: {e}')

    @property
    def sources(self) -> list[CalibratorSource]:
        """List of all `CalibratorSource` objects currently in the catalog."""
        return self._sources

    @property
    def n_sources(self) -> int:
        """Number of sources currently in the catalog."""
        return len(self._sources)

    def get_source(self, name: str) -> Optional[CalibratorSource]:
        """Get a source by name, IVS name, or other names, all case-insensitive.

        Tries, in order: exact match on primary name, exact match on IVS name, exact match on
        any of `source.other_names`, then finally a partial (substring, either direction)
        match on primary name and on IVS name.

        Inputs
        - name : str
            Name, IVS name, alias, or substring thereof, to search for.

        Returns
        - CalibratorSource or None
            The first matching source, or None if no match is found.
        """
        name_upper = name.upper()

        # First try exact matches (case-insensitive)
        for source_name, source in self._name_index.items():
            if source_name.upper() == name_upper:
                return source

        for ivs_name, source in self._ivsname_index.items():
            if ivs_name.upper() == name_upper:
                return source

        # Then search in other names
        for source in self._sources:
            if source.other_names:
                for other_name in source.other_names:
                    if other_name.upper() == name_upper:
                        return source

        # Finally, try partial matches (e.g., if name is contained in a longer name)
        for source_name, source in self._name_index.items():
            if name_upper in source_name.upper() or source_name.upper() in name_upper:
                return source

        for ivs_name, source in self._ivsname_index.items():
            if name_upper in ivs_name.upper() or ivs_name.upper() in name_upper:
                return source

        return None

    def calibrators_only(self) -> Self:
        """Returns a new RFCCatalog containing only sources flagged as calibrators.

        Returns
        - RFCCatalog
            A new catalog instance (same `min_flux`/`band`/`catalog_filename` settings) whose
            sources are the subset of `self.sources` with `is_calibrator == True`.
        """
        new_catalog = object.__new__(self.__class__)
        new_catalog._sources = [s for s in self._sources if s.is_calibrator]
        new_catalog._min_flux = self._min_flux
        new_catalog._band = self._band
        new_catalog._catalog_filename = self._catalog_filename
        new_catalog._name_index = {s.name: s for s in new_catalog._sources}
        new_catalog._ivsname_index = {s.ivsname: s for s in new_catalog._sources}
        if new_catalog._sources:
            new_catalog._ra_arr = np.array([s.ra_deg for s in new_catalog._sources], dtype=np.float64)
            new_catalog._dec_arr = np.array([s.dec_deg for s in new_catalog._sources], dtype=np.float64)
        else:
            new_catalog._ra_arr = np.array([])
            new_catalog._dec_arr = np.array([])
        return new_catalog

    def brighter_than(self, flux: float, band: Optional[str] = None) -> Self:
        """Returns a new RFCCatalog containing only sources above a given unresolved flux.

        Inputs
        - flux : float
            Minimum unresolved flux density (Jy) a source must have at `band` to be kept.
        - band : str [OPTIONAL]
            RFC band letter code to check the flux at. If not given, uses the catalog's own
            `band` (the one used at load time).

        Returns
        - RFCCatalog
            A new catalog instance (same `min_flux`/`band`/`catalog_filename` settings) whose
            sources are the subset of `self.sources` with `unresolved_flux(band) >= flux`.
        """
        check_band = band if band is not None else self._band
        new_catalog = object.__new__(self.__class__)
        new_catalog._sources = [s for s in self._sources if s.unresolved_flux(check_band) >= flux]
        new_catalog._min_flux = self._min_flux
        new_catalog._band = self._band
        new_catalog._catalog_filename = self._catalog_filename
        new_catalog._name_index = {s.name: s for s in new_catalog._sources}
        new_catalog._ivsname_index = {s.ivsname: s for s in new_catalog._sources}
        if new_catalog._sources:
            new_catalog._ra_arr = np.array([s.ra_deg for s in new_catalog._sources], dtype=np.float64)
            new_catalog._dec_arr = np.array([s.dec_deg for s in new_catalog._sources], dtype=np.float64)
        else:
            new_catalog._ra_arr = np.array([])
            new_catalog._dec_arr = np.array([])
        return new_catalog

    def _get_coord_arrays(self) -> tuple[np.ndarray, np.ndarray]:
        """Returns the cached RA/Dec coordinate arrays for all sources in the catalog.

        Returns
        - tuple[numpy.ndarray, numpy.ndarray]
            (ra_deg_array, dec_deg_array), in the same order as `self.sources`.
        """
        return self._ra_arr, self._dec_arr


def _angular_separation(ra1_deg: float, dec1_deg: float, ra2_arr: np.ndarray, dec2_arr: np.ndarray) -> np.ndarray:
    """Vectorized angular separation (haversine formula) between one point and an array of points.

    Inputs
    - ra1_deg : float
        Right ascension of the reference point, in degrees.
    - dec1_deg : float
        Declination of the reference point, in degrees.
    - ra2_arr : numpy.ndarray
        Right ascensions of the target points, in degrees.
    - dec2_arr : numpy.ndarray
        Declinations of the target points, in degrees. Must be the same shape as `ra2_arr`.

    Returns
    - numpy.ndarray
        Angular separation in degrees between the reference point and each target point,
        same shape as `ra2_arr`.
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
    once per time step, then transforming all sources vectorized.

    Inputs
    - ra_rad : numpy.ndarray
        Right ascensions of the sources, in radians, shape (n_sources,).
    - dec_rad : numpy.ndarray
        Declinations of the sources, in radians, shape (n_sources,).
    - times : astropy.time.Time
        Time array (length n_times) at which to evaluate the positions.
    - station : Station
        Station providing the observing location (`station.location`).

    Returns
    - tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]
        (elevation_deg, azimuth_deg, ha_hours), each shaped (n_times, n_sources). `ha_hours`
        is wrapped to [0, 24).
    """
    loc = station.location
    lon_rad, lat_rad = loc.lon.rad, loc.lat.rad
    height_m = loc.height.to(u.m).value
    utc1, utc2 = times.utc.jd1, times.utc.jd2
    dut1 = times.delta_ut1_utc
    n_times = len(times)
    n_src = len(ra_rad)
    elev_out = np.empty((n_times, n_src))
    az_out = np.empty((n_times, n_src))
    ha_out = np.empty((n_times, n_src))
    for t in range(n_times):
        astrom, eo = erfa.apco13(utc1[t], utc2[t], dut1[t], lon_rad, lat_rad, height_m, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        ri, di = erfa.atciq(ra_rad, dec_rad, 0.0, 0.0, 0.0, 0.0, astrom)
        az, zen, ha, _dec, _ra = erfa.atioq(ri, di, astrom)
        elev_out[t] = np.degrees(np.pi / 2.0 - zen)
        az_out[t] = np.degrees(az)
        ha_out[t] = (np.degrees(ha) / 15.0) % 24.0
    return elev_out, az_out, ha_out


def _station_observable_mask(elev: np.ndarray, az: np.ndarray, ha_hours: np.ndarray,
                             dec_deg: np.ndarray, station) -> np.ndarray:
    """Build boolean observable mask (n_times, n_sources) respecting station mount constraints.

    Uses the same constraint logic as astroplan: for ALTAZ mounts checks azimuth and elevation
    limits; for EQUAT mounts checks hour angle (in [0,24h) with wrapping), declination, and
    a minimum 5-degree elevation. If the station defines a local horizon, also requires the
    elevation to be above the interpolated horizon at each azimuth.

    Inputs
    - elev : numpy.ndarray
        Elevation in degrees, shape (n_times, n_sources).
    - az : numpy.ndarray
        Azimuth in degrees, same shape as `elev`.
    - ha_hours : numpy.ndarray
        Hour angle in hours (wrapped to [0, 24)), same shape as `elev`.
    - dec_deg : numpy.ndarray
        Declination in degrees, shape (n_sources,).
    - station : Station
        Station providing the mount constraints (`station.mount`) and optional local horizon
        (`station.horizon`, `station.horizon_min_elevation`).

    Returns
    - numpy.ndarray
        Boolean mask, shape (n_times, n_sources), True where the source is observable.
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

    Inputs
    - stations : Stations
        Network of stations to check visibility for.
    - times : astropy.time.Time
        Time array sampling the observation window.
    - min_elevation : astropy.units.Quantity [OPTIONAL]
        Minimum elevation a source must reach to be considered visible. Defaults to 20 deg.
    - min_flux : astropy.units.Quantity [OPTIONAL]
        Minimum unresolved flux density (in Jy) used to build the default catalog, when
        `catalog` is not given. Defaults to 1.0 Jy. Ignored if `catalog` is provided.
    - catalog : RFCCatalog [OPTIONAL]
        Pre-loaded catalog to search. If not given, a new `RFCCatalog(min_flux=min_flux,
        band='c')` is created.
    - require_all_stations : bool [OPTIONAL]
        If True (default), only sources visible by ALL stations for the entire time range are
        returned. If False, sources visible by at least one station at any time are returned,
        together with per-source antenna-visibility details.

    Returns
    - tuple[list[CalibratorSource], list[float], Optional[list[tuple[int, int, bool, float]]]]
        (sources, min_elevations_deg, antenna_visibility). `sources` and `min_elevations_deg`
        are sorted best-first (highest minimum elevation, or best visibility category when
        `require_all_stations` is False). `antenna_visibility` is None when
        `require_all_stations` is True; otherwise a list of
        (n_stations_visible, n_stations_total, visible_all_times, min_elevation_deg_across_
        visible_stations) tuples, one per returned source.
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

        if len(visible_idx) == 0:
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

        if len(visible_idx) == 0:
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
    """Find calibrator sources near a target source, sorted by increasing separation.

    Inputs
    - source : CalibratorSource or Source
        Target source to search around.
    - max_separation : astropy.units.Quantity [OPTIONAL]
        Maximum angular separation from `source` to include. Defaults to 5 deg.
    - catalog : RFCCatalog [OPTIONAL]
        Pre-loaded catalog to search. If not given, a new default `RFCCatalog()` is created.
    - n_sources : int [OPTIONAL]
        Maximum number of sources to return. If not given, all matches are returned.

    Returns
    - list[tuple[CalibratorSource, float]]
        List of (source, separation_deg) tuples, sorted by increasing separation, excluding
        the target itself (separation == 0 is filtered out) and truncated to `n_sources` if
        given.
    """
    if catalog is None:
        catalog = RFCCatalog()

    max_sep_deg = max_separation.to(u.deg).value if hasattr(max_separation, 'to') else max_separation
    ra_arr, dec_arr = catalog._get_coord_arrays()
    if len(ra_arr) == 0:
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

    Prioritises the brightest, most compact unresolved emission among the closest sources.
    The scoring formula is ``score = unresolved_flux * compactness / (1 + separation_deg)``,
    where ``compactness = min(unresolved_flux / resolved_flux, 1)``, so that a bright, nearby,
    compact source is preferred.

    Inputs
    - target : Source
        The target source to find a phase calibrator for.
    - band : str
        Observing band in RFC letter code (s/c/x/u/k) or wavelength string (e.g. '6cm').
    - max_separation : astropy.units.Quantity [OPTIONAL]
        Maximum angular separation from target. Defaults to 5 deg.
    - catalog : RFCCatalog [OPTIONAL]
        Pre-loaded catalog. If None, a new one is created with min_flux=0 at the given band.

    Returns
    - CalibratorSource or None
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

    Inputs
    - target : Source
        The target source.
    - phase_cal : Source
        The already-selected phase calibrator.
    - band : str
        Observing band (RFC letter or wavelength string).
    - max_separation : astropy.units.Quantity [OPTIONAL]
        Maximum angular separation from target. Defaults to 5 deg.
    - catalog : RFCCatalog [OPTIONAL]
        Pre-loaded catalog. If None, a new one is created with min_flux=0 at the given band.

    Returns
    - CalibratorSource or None
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

    Also accepts bare RFC band letters ('s', 'c', 'x', 'u', 'k') unchanged, and unknown
    wavelength strings (rounded to the nearest known one via `_round_to_nearest_wavelength`).
    Falls back to 'c' if the band cannot be matched at all.

    Inputs
    - band : str
        Wavelength string (e.g. '6cm', '18cm') or RFC letter.

    Returns
    - str
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


def main_fringe():
    """CLI entry point ('planobs_fringefinder') to find and print fringe finder candidates
    for a VLBI observation, given a network/station selection and a time range.

    Parses command-line arguments, resolves the requested stations/networks, searches the RFC
    catalog via `get_fringe_finder_sources`, and prints the results as a rich table or as JSON
    (`--json`). Exits the process (`sys.exit`) on both success and error paths; does not
    return a value.
    """
    usage = "%(prog)s [-h] OPTIONS"
    description = "Find fringe finder sources for VLBI observations"
    parser = argparse.ArgumentParser(description=description, prog="planobs_fringefinder", usage=usage,
                                  formatter_class=RawTextRichHelpFormatter)
    parser.add_argument('-n', '--network', type=str, nargs='+',
                        help="The VLBI network(s) that will participate in the observation. "
                        "It will take the default stations in each network.")
    parser.add_argument('-s', '--stations', type=str, nargs='+',
                        help="List of antenna codenames or names that will participate in the observation.")
    parser.add_argument('-t', '--starttime', type=str, required=True,
                        help="Start of the observation in format 'YYYY-MM-DD HH:MM' (UTC).")
    parser.add_argument('-d', '--duration', type=float, required=True, help="Duration of the observation in hours.")
    parser.add_argument('--min-flux', type=float, default=0.5,
                        help="Minimum unresolved flux threshold in Jy (default: 0.5).")
    parser.add_argument('--min-elevation', type=float, default=20.0,
                        help="Minimum elevation in degrees (default: 20).")
    parser.add_argument('-l', '--max-lines', type=int, default=20,
                        help="Maximum number of sources to return (default: 20).")
    parser.add_argument('--require-all', action='store_true', default=False,
                        help="Require source to be visible by ALL stations (default: False).")
    parser.add_argument('-b', '--band', type=str, default=None,
                        help="Observing band for flux display (e.g., '18cm', '6cm'). "
                             "If not provided, shows flux for all available bands.")
    parser.add_argument('--station-catalog', type=str, default=None, help="Path to custom station catalog file.")
    parser.add_argument('--json', action='store_true', default=False,
                        help="Output results in JSON format instead of a table.")
    args = parser.parse_args()
    if args.network is None and args.stations is None:
        error_msg = "You need to provide at least a VLBI network or a list of antennas."
        if args.json:
            print(json.dumps({"error": error_msg}, indent=2))
        else:
            rprint(f"[bold red]{error_msg}[/bold red]")
        sys.exit(1)

    obs._STATIONS = obs.Stations(filename=args.station_catalog)
    stations_list = []
    if args.network is not None:
        try:
            for n in args.network:
                network = obs._NETWORKS[n]
                for s in network.station_codenames:
                    if s not in stations_list:
                        stations_list.append(s)
        except KeyError:
            unknown_networks = [n for n in args.network if n not in obs._NETWORKS]
            n_networks = len(unknown_networks)
            error_msg = (f"The network{'s' if n_networks > 1 else ''} {', '.join(unknown_networks)} "
                         f"{'are' if n_networks > 1 else 'is'} not known.")
            if args.json:
                print(json.dumps({"error": error_msg}, indent=2))
            else:
                rprint(f"[bold red]{error_msg}[/bold red]")
            sys.exit(1)

    for s in args.stations or []:
        try:
            # Try case-sensitive lookup first
            a_station = obs._STATIONS[s.strip()].codename
        except KeyError:
            try:
                # Try case-insensitive lookup by searching through all stations
                s_upper = s.strip().upper()
                found_station = None
                for codename in obs._STATIONS.station_codenames:
                    if codename.upper() == s_upper:
                        found_station = obs._STATIONS[codename].codename
                        break

                if found_station is None:
                    # Also try full station names (case insensitive)
                    for name in obs._STATIONS.station_names:
                        if name.upper() == s_upper:
                            found_station = obs._STATIONS[name].codename
                            break

                if found_station is None:
                    raise KeyError(f"Station {s} not found")

                a_station = found_station
            except KeyError:
                error_msg = f"The station {s} is not known."
                if args.json:
                    print(json.dumps({"error": error_msg}, indent=2))
                else:
                    rprint(f"[bold red]{error_msg}[/bold red]")
                sys.exit(1)

        if a_station not in stations_list:
            stations_list.append(a_station)

    stations_obj = obs._STATIONS.filter_antennas(stations_list)
    if not stations_obj:
        error_msg = "No valid antennas have been selected."
        if args.json:
            print(json.dumps({"error": error_msg}, indent=2))
        else:
            rprint(f"[bold red]{error_msg}[/bold red]")
        sys.exit(1)

    times = Time(args.starttime, scale='utc') + np.arange(0, args.duration + 0.1, 0.1) * u.hour
    sources, min_elevs, antenna_visibility = get_fringe_finder_sources(stations_obj, times,
                                               min_elevation=args.min_elevation * u.deg,
                                               min_flux=args.min_flux * u.Jy,
                                               require_all_stations=args.require_all)
    if not sources:
        err_msg = (f"No fringe finder candidates found above {args.min_elevation} degrees elevation "
                   f"and with a unresolved flux above {args.min_flux} Jy.")
        result = {"error": err_msg}
        if args.json:
            print(json.dumps(result, indent=2))
        else:
            rprint(f"[bold red]{err_msg}[/bold red]")
        sys.exit(0)

    max_display = min(args.max_lines, len(sources))
    sources_slice = sources[:max_display]
    min_elevs_slice = min_elevs[:max_display]

    result_data = []
    if antenna_visibility is not None:
        for i in range(max_display):
            src = sources_slice[i]
            min_elev = min_elevs_slice[i]
            visible_count, total_count, visible_all_times, min_elev_all = antenna_visibility[i]

            if visible_count == total_count and visible_all_times:
                visibility_text = "all ant. all time"
            elif visible_count == total_count:
                visibility_text = "all ant. partial time"
            else:
                visibility_text = f"{visible_count}/{total_count} ant."

            # Get flux information for the specified band
            if args.band:
                total_flux, unresolved_flux = src.get_flux_at_band(args.band)
            else:
                # Use first available band for total flux display
                total_flux = float(np.max(src.flux_resolved)) if np.any(src.flux_resolved > 0) else 0.0
                unresolved_flux = float(np.max(src.flux_unresolved)) if np.any(src.flux_unresolved > 0) else 0.0

            result_data.append({"name": src.name, "ivs_name": src.ivsname,
                            "min_elevation_deg": min_elev if min_elev > 0.0 else 0,
                            "total_flux_jy": total_flux, "unresolved_flux_jy": unresolved_flux,
                            "bands": src.get_observed_bands(),
                            "astrogeo_url": src.get_astrogeo_link(),
                            "antenna_visibility": visibility_text})
    else:
        for i in range(max_display):
            src = sources_slice[i]
            min_elev = min_elevs_slice[i]

            # Get flux information for the specified band
            if args.band:
                total_flux, unresolved_flux = src.get_flux_at_band(args.band)
            else:
                # Use first available band for total flux display
                total_flux = float(np.max(src.flux_resolved)) if np.any(src.flux_resolved > 0) else 0.0
                unresolved_flux = float(np.max(src.flux_unresolved)) if np.any(src.flux_unresolved > 0) else 0.0

            result_data.append({"name": src.name, "ivs_name": src.ivsname,
                            "min_elevation_deg": min_elev if min_elev > 0.0 else 0,
                            "total_flux_jy": total_flux, "unresolved_flux_jy": unresolved_flux,
                            "bands": src.get_observed_bands(),
                            "astrogeo_url": src.get_astrogeo_link()})

    result = {"min_elevation_deg": args.min_elevation, "min_flux_jy": args.min_flux,
              "require_all_stations": args.require_all, "sources": result_data,
              "total_found": len(sources), "shown": max_display}

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        rprint(f"\n[bold green]Found {len(sources)} fringe finder candidates above "
               f"{args.min_elevation} degrees elevation and with a unresolved flux "
               f"above {args.min_flux} Jy:[/bold green]")

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
            src = sources_slice[i]
            min_elev = min_elevs_slice[i]

            # Get flux information for the specified band
            if args.band:
                total_flux, unresolved_flux = src.get_flux_at_band(args.band)
            else:
                # Use first available band for total flux display
                total_flux = float(np.max(src.flux_resolved)) if np.any(src.flux_resolved > 0) else 0.0
                unresolved_flux = float(np.max(src.flux_unresolved)) if np.any(src.flux_unresolved > 0) else 0.0

            row = [src.name, src.ivsname, f"{min_elev if min_elev > 0.0 else 0:>6.1f}",
                   f"{total_flux:>8.2f}" if total_flux > 0 else "N/A",
                   f"{unresolved_flux:>8.2f}" if unresolved_flux > 0 else "N/A",
                   src.get_observed_bands(), f"[link={src.get_astrogeo_link()}]AstroGeo[/link]"]

            if antenna_visibility is not None:
                visible_count, total_count, visible_all_times, min_elev_all = antenna_visibility[i]
                if visible_count == total_count and visible_all_times:
                    visibility_text = "all, all time"
                elif visible_count == total_count:
                    visibility_text = "all, partial time"
                else:
                    visibility_text = f"{visible_count}/{total_count} antennas"
                row.append(visibility_text)

            table.add_row(*row)

        rprint(table)
        if len(sources) > max_display:
            rprint(f"\n... and {len(sources) - max_display} more sources.")
    sys.exit(0)


def _target_from_personal_catalog(catalog_file: str, name: str) -> Optional[Source]:
    """Looks up `name` in a personal source catalog (toml file, as used by 'planobs observe -sc').

    Inputs
    - catalog_file : str
        Path to the personal source catalog (toml).
    - name : str
        Block name or source name to look up.

    Returns
    - Source or None
        The (first) target source of the matching block, the source with the given name, or
        None if not found.

    Raises
    - FileNotFoundError: If catalog_file does not exist.
    """
    catalog = SourceCatalog(catalog_file)
    if name in catalog.blocknames:
        block = catalog[name]
        block_targets = block.sources(SourceType.TARGET) or block.sources()
        return block_targets[0] if block_targets else None

    return catalog.sources(include_calibrators=True).get(name)


def main_phasecal():
    """CLI entry point ('planobs_phasecal') to find and print phase calibrator candidates
    near a target source.

    Parses command-line arguments, resolves the target source (optionally via a personal
    source catalog first, then the RFC catalog, then a free-form parse via
    `Source.source_from_str`), searches nearby sources via `get_nearby_sources`, and prints
    the results as a rich table or as JSON (`--json`). Exits the process (`sys.exit`) on both
    success and error paths; does not return a value.
    """
    usage = "%(prog)s [-h] OPTIONS"
    description = "Find phase calibrator sources near a target source"
    parser = argparse.ArgumentParser(description=description, prog="planobs_phasecal", usage=usage,
                                  formatter_class=RawTextRichHelpFormatter)
    parser.add_argument('-t', '--target', type=str, required=True,
                        help="Target source name (J2000 or IVS name from RFC catalog).")
    parser.add_argument('--max-separation', type=float, default=5.0,
                        help="Maximum angular separation in degrees (default: 5.0).")
    parser.add_argument('--min-flux', type=float, default=0.1,
                        help="Minimum unresolved flux threshold in Jy (default: 0.1).")
    parser.add_argument('-n', '--n-sources', type=int, default=None,
                        help="Maximum number of sources to return (default: all).")
    parser.add_argument('-b', '--band', type=str, default=None,
                        help="Observing band for flux display (e.g., '18cm', '6cm'). "
                             "If not provided, shows flux for all available bands.")
    parser.add_argument('--catalog-file', type=str, default=None,
                        help="Path to custom RFC catalog file.")
    parser.add_argument('-sc', '--source-catalog', '--sc', type=str, default=None,
                        help="Input file containing the personal source catalog (toml).\n"
                        "If provided, '--target' is first looked up there (block or source name).")
    parser.add_argument('--json', action='store_true', default=False,
                        help="Output results in JSON format instead of a table.")
    args = parser.parse_args()

    catalog = RFCCatalog(catalog_filename=args.catalog_file, band='c', min_flux=args.min_flux, include_missing=True)
    target = None
    if args.source_catalog is not None:
        try:
            target = _target_from_personal_catalog(args.source_catalog, args.target)
        except FileNotFoundError:
            error_msg = f"Source catalog file not found: {args.source_catalog}"
            if args.json:
                print(json.dumps({"error": error_msg}, indent=2))
            else:
                rprint(f"[bold red]{error_msg}[/bold red]")
            sys.exit(1)
        if target is None:
            rprint(f"[yellow]'{args.target}' not found in {args.source_catalog} — "
                   "falling back to the RFC catalog/online lookup.[/yellow]")

    if not target:
        target = catalog.get_source(args.target)
    if not target:
        try:
            target = Source.source_from_str(args.target, source_type=SourceType.TARGET)
        except Exception:
            error_msg = f"Target source '{args.target}' could not be parsed or found in the catalogs."
            if args.json:
                print(json.dumps({"error": error_msg}, indent=2))
            else:
                rprint(f"[bold red]{error_msg}[/bold red]")
            sys.exit(1)

    nearby = get_nearby_sources(target, max_separation=args.max_separation * u.deg,
                               catalog=catalog, n_sources=args.n_sources)
    if not nearby:
        result = {"error": f"No phase calibrator candidates found near {target.name} "
                           f"({target.coord.to_string('hmsdms')})."}
        if args.json:
            print(json.dumps(result, indent=2))
        else:
            rprint(f"[bold red]No phase calibrator candidates found near {target.name} "
                   f"({target.coord.to_string('hmsdms')}).[/bold red]")
        sys.exit(1)

    result_data = []
    for src, sep in nearby:
        # Get flux information for the specified band
        if args.band:
            total_flux, unresolved_flux = src.get_flux_at_band(args.band)
        else:
            # Use first available band for total flux display
            total_flux = float(np.max(src.flux_resolved)) if np.any(src.flux_resolved > 0) else 0.0
            unresolved_flux = float(np.max(src.flux_unresolved)) if np.any(src.flux_unresolved > 0) else 0.0

        result_data.append({"name": src.name, "ivs_name": src.ivsname, "separation_deg": sep,
                            "total_flux_jy": total_flux, "unresolved_flux_jy": unresolved_flux,
                            "bands": src.get_observed_bands(), "astrogeo_url": src.get_astrogeo_link()})

    result = {"target_name": target.name, "target_coordinates": target.coord.to_string('hmsdms'),
              "max_separation_deg": args.max_separation, "min_flux_jy": args.min_flux,
              "sources": result_data, "total_found": len(nearby)}

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        rprint(f"\n[bold green]Found {len(nearby)} phase calibrator candidates near {target.name} "
               f"({target.coord.to_string('hmsdms')}):[/bold green]")

        table = Table(show_header=True, header_style="bold", show_lines=False, box=box.SIMPLE)
        table.add_column("Name", style="", width=17)
        table.add_column("IVS Name", style="", width=10)
        table.add_column("Separation (deg)", justify="right", style="", width=12)
        table.add_column("Total flux (Jy)", justify="right", style="", width=12)
        table.add_column("Unresolved (Jy)", justify="right", style="", width=13)
        table.add_column("Bands", justify="right", style="", width=10)
        table.add_column("url", style="", width=10)

        for src, sep in nearby:
            # Get flux information for the specified band
            if args.band:
                total_flux, unresolved_flux = src.get_flux_at_band(args.band)
            else:
                # Use first available band for total flux display
                total_flux = float(np.max(src.flux_resolved)) if np.any(src.flux_resolved > 0) else 0.0
                unresolved_flux = float(np.max(src.flux_unresolved)) if np.any(src.flux_unresolved > 0) else 0.0

            table.add_row(src.name, src.ivsname, f"{sep:.2f}",
                          f"{total_flux:>8.2f}" if total_flux > 0 else "N/A",
                          f"{unresolved_flux:>8.2f}" if unresolved_flux > 0 else "N/A",
                          src.get_observed_bands(),
                          f"[link={src.get_astrogeo_link()}]AstroGeo[/link]")

        rprint(table)
    sys.exit(0)
