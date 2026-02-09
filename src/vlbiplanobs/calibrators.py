"""Calibrator source module for finding fringe finders and nearby sources."""

import sys
import argparse
from typing import Optional, Self
from importlib import resources
import numpy as np
from astropy import units as u, coordinates as coord
from astropy.time import Time
from astroplan import FixedTarget
from rich import print as rprint, box
from rich.table import Table
from rich_argparse import RawTextRichHelpFormatter
from urllib import parse

from .sources import Source, SourceType
from .stations import Stations
from . import observation as obs

_RFC_BANDS = {'l': 's', 's': 's', 'c': 'c', 'm': 'c', 'x': 'x', 'u': 'u', 'k': 'k', 'q': 'k'}
_DEFAULT_MIN_ELEVATION = 20 * u.deg
_DEFAULT_MIN_FLUX = 1.0 * u.Jy
_BAND_INDEX = {'s': 0, 'c': 1, 'x': 2, 'u': 3, 'k': 4}
_WAVELENGTH_BANDS = {'18cm': 's', '21cm': 's', '13cm': 'c', '6cm': 'c', '5cm': 'c', 
                       '3.6cm': 'x', '2cm': 'u', '1.3cm': 'k', '0.7cm': 'k'}


def _round_to_nearest_wavelength(band_str: str) -> str:
    """Round an unknown wavelength to the nearest known wavelength."""
    if band_str in _WAVELENGTH_BANDS:
        return band_str
    if not band_str.lower().endswith('cm'):
        return band_str
    return min(_WAVELENGTH_BANDS.keys(), key=lambda x: abs(float(x[:-2]) - float(band_str[:-2])))


class CalibratorSource(Source):
    """Represents a calibrator source from the RFC catalog."""
    __slots__ = ('ivsname', 'n_observations', 'flux_resolved', 'flux_unresolved', 'is_calibrator')

    def __init__(self, name: str, ivsname: str, ra_deg: float, dec_deg: float, n_observations: int, 
                 flux_resolved: np.ndarray, flux_unresolved: np.ndarray, is_calibrator: bool):
        super().__init__(name=name, coordinates=coord.SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg), 
                     source_type=SourceType.PHASECAL)
        self.ivsname = ivsname
        self.n_observations = n_observations
        self.flux_resolved = flux_resolved
        self.flux_unresolved = flux_unresolved
        self.is_calibrator = is_calibrator

    @property
    def ra_deg(self) -> float:
        return float(self.coord.ra.deg)

    @property
    def dec_deg(self) -> float:
        return float(self.coord.dec.deg)

    def unresolved_flux(self, band: str) -> float:
        idx = _BAND_INDEX.get(band)
        return float(self.flux_unresolved[idx]) if idx is not None else 0.0

    def resolved_flux(self, band: str) -> float:
        idx = _BAND_INDEX.get(band)
        return float(self.flux_resolved[idx]) if idx is not None else 0.0

    def get_flux_at_band(self, band: str) -> tuple[float, float]:
        """Get resolved and unresolved flux at a specific band with interpolation."""
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
        """Interpolate flux for a band using nearby bands."""
        # Find the wavelength for the target band
        target_wavelength = None
        for wl, band in _WAVELENGTH_BANDS.items():
            if band == target_band:
                target_wavelength = float(wl[:-2])
                break
        
        if target_wavelength is None:
            return 0.0, 0.0
        
        # Vectorized approach: collect all available data at once
        available_data = []
        for band, idx in _BAND_INDEX.items():
            if idx < len(self.flux_unresolved) and self.flux_unresolved[idx] > 0:
                for wl, wl_band in _WAVELENGTH_BANDS.items():
                    if wl_band == band:
                        available_data.append((float(wl[:-2]), float(self.flux_resolved[idx]), float(self.flux_unresolved[idx])))
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
            if res1 > 0 and res2 > 0:
                log_res1, log_res2 = np.log10(res1), np.log10(res2)
                log_wl1, log_wl2 = np.log10(wl1), np.log10(wl2)
                log_target_wl = np.log10(target_wavelength)
                log_res_interp = log_res1 + (log_res2 - log_res1) * (log_target_wl - \
                    log_wl1) / (log_wl2 - log_wl1)
                res_interp = 10**log_res_interp
            else:
                res_interp = (res1 + res2) / 2
            
            if unres1 > 0 and unres2 > 0:
                log_unres1, log_unres2 = np.log10(unres1), np.log10(unres2)
                log_unres_interp = log_unres1 + (log_unres2 - log_unres1) * (log_target_wl - log_wl1) / (log_wl2 - log_wl1)
                unres_interp = 10**log_unres_interp
            else:
                unres_interp = (unres1 + unres2) / 2
        else:
            res_interp = (res1 + res2) / 2
            unres_interp = (unres1 + unres2) / 2
        
        return res_interp, unres_interp

    def get_skycoord(self) -> coord.SkyCoord:
        return self.coord

    def get_astrogeo_link(self) -> str:
        source_coord_str = parse.quote("ra={:02.0f}:{:02.0f}:{:06.3f}&dec={:+03.0f}:{:02.0f}:{:06.3f}&num_sou=1&format=html".format(
            *self.coord.ra.hms, *self.coord.dec.dms), safe='=&')
        return f"http://astrogeo.org/cgi-bin/calib_search_form.csh?{source_coord_str}"

    def get_observed_bands(self) -> str:
        bands = []
        for band in ['s', 'c', 'x', 'u', 'k']:
            idx = _BAND_INDEX[band]
            if self.flux_unresolved[idx] > 0 or self.flux_resolved[idx] > 0:
                bands.append(band.upper())
        return ','.join(bands)


class RFCCatalog:
    """RFC (Radio Fundamental Catalog) of VLBI calibrator sources."""
    __slots__ = ('_sources', '_min_flux', '_band', '_catalog_filename', 
                 '_name_index', '_ivsname_index', '_ra_arr', '_dec_arr')

    def __init__(self, catalog_filename: Optional[str] = None,
                min_flux: u.Quantity = _DEFAULT_MIN_FLUX, band: str = 'c'):
        self._sources: list[CalibratorSource] = []
        self._min_flux = min_flux.to(u.Jy).value if hasattr(min_flux, 'to') else min_flux
        self._band = band
        self._catalog_filename = catalog_filename
        self._name_index: dict[str, CalibratorSource] = {}
        self._ivsname_index: dict[str, CalibratorSource] = {}
        self._ra_arr: Optional[np.ndarray] = None
        self._dec_arr: Optional[np.ndarray] = None
        self._load_catalog()

    def _get_catalog_path(self) -> str:
        if self._catalog_filename is not None:
            return self._catalog_filename
        rfc_files = tuple(r.name for r in resources.files('vlbiplanobs.data').iterdir() 
                       if r.is_file() and 'rfc' in r.name and r.name.endswith('.txt'))
        if len(rfc_files) == 0:
            raise FileNotFoundError('No RFC catalog files found in the data directory.')
        with resources.as_file(resources.files('vlbiplanobs.data').joinpath(sorted(rfc_files)[-1])) as rfcfile:
            return str(rfcfile)

    def _load_catalog(self) -> None:
        catalog_path = self._get_catalog_path()
        try:
            band_idx = _BAND_INDEX[self._band]
            with open(catalog_path, 'rt') as fin:
                lines = fin.readlines()
            
            valid_lines = [line for line in lines 
                          if line and line[0] not in '#NU' and len(line) >= 100]
            
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
                    if flux_array[band_idx, 1] < self._min_flux:
                        continue
                except (ValueError, IndexError):
                    continue
                
                try:
                    ra_deg = (float(cols[3]) + float(cols[4]) / 60.0 + float(cols[5]) / 3600.0) * 15.0
                    dec_deg = (1.0 if cols[6][0] != '-' else -1.0) * (abs(float(cols[6])) + float(cols[7]) / 60.0 + float(cols[8]) / 3600.0)
                except (ValueError, IndexError):
                    continue
                
                parsed_data.append({'name': cols[2], 'ivsname': cols[1], 'ra': ra_deg, 'dec': dec_deg,
                    'n_obs': int(cols[12]) if cols[12].lstrip('-').isdigit() else 0,
                    'flux_resolved': flux_array[:, 0], 'flux_unresolved': flux_array[:, 1],
                    'is_cal': cols[0] == 'C'})
            
            self._sources = [CalibratorSource(d['name'], d['ivsname'], d['ra'], d['dec'], d['n_obs'],
                d['flux_resolved'], d['flux_unresolved'], d['is_cal']) for d in parsed_data]
            
            self._name_index = {s.name: s for s in self._sources}
            self._ivsname_index = {s.ivsname: s for s in self._sources}
            
            if self._sources:
                coords = np.array([(s.ra_deg, s.dec_deg) for s in self._sources], dtype=np.float64)
                self._ra_arr = coords[:, 0].copy()
                self._dec_arr = coords[:, 1].copy()
        except FileNotFoundError:
            raise FileNotFoundError(f'RFC catalog file not found: {catalog_path}')
        except Exception as e:
            raise RuntimeError(f'Error parsing RFC catalog file: {e}')

    @property
    def sources(self) -> list[CalibratorSource]:
        return self._sources

    @property
    def n_sources(self) -> int:
        return len(self._sources)

    def get_source(self, name: str) -> Optional[CalibratorSource]:
        return self._name_index.get(name) or self._ivsname_index.get(name)

    def calibrators_only(self) -> Self:
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
        return self._ra_arr, self._dec_arr


def _angular_separation(ra1_deg: float, dec1_deg: float, ra2_arr: np.ndarray, dec2_arr: np.ndarray) -> np.ndarray:
    """Vectorized angular separation calculation."""
    ra1_rad, dec1_rad = np.radians(ra1_deg), np.radians(dec1_deg)
    ra2_rad, dec2_rad = np.radians(ra2_arr), np.radians(dec2_arr)
    
    # Vectorized haversine formula
    dlon = ra2_rad - ra1_rad
    dlat = dec2_rad - dec1_rad
    
    a = np.sin(dlat / 2.0) ** 2 + np.cos(dec1_rad) * np.cos(dec2_rad) * np.sin(dlon / 2.0) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(a)))


def get_fringe_finder_sources(stations: Stations, times: Time, band: str, 
                           min_elevation: u.Quantity = _DEFAULT_MIN_ELEVATION, 
                           min_flux: u.Quantity = _DEFAULT_MIN_FLUX, 
                           catalog: Optional[RFCCatalog] = None, 
                           require_all_stations: bool = True) -> tuple[list[CalibratorSource], list[float]]:
    """Find fringe finder sources visible by the given stations.
    
    Args:
        stations: List of observing stations
        times: Observation times
        band: Observing band (RFC band, wavelength, or None for all bands)
        min_elevation: Minimum elevation requirement
        min_flux: Minimum flux requirement
        catalog: RFC catalog (created if None)
        require_all_stations: If True, source must be visible by all stations
        
    Returns:
        Tuple of (sources, minimum_elevations)
    """
    # Convert wavelength to RFC band if needed
    if band in _WAVELENGTH_BANDS:
        band = _WAVELENGTH_BANDS[band]
    rfc_band = _RFC_BANDS.get(band.lower(), band.lower())
    bands_to_check = [rfc_band]
    
    if catalog is None:
        catalog = RFCCatalog(min_flux=min_flux, band=bands_to_check[0])
    
    if len(stations.stations) == 0:
        return [], []
    
    min_el_deg = min_elevation.to(u.deg).value if hasattr(min_elevation, 'to') else min_elevation
    min_flux_val = min_flux.value if hasattr(min_flux, 'value') else min_flux
    
    # Pre-filter sources by flux to reduce computation
    valid_sources = []
    for src in catalog.sources:
        _, unresolved_flux = src.get_flux_at_band(bands_to_check[0])
        if unresolved_flux >= min_flux_val:
            valid_sources.append(src)
    
    if not valid_sources:
        return [], []
    
    visible_sources = []
    min_elevations = []
    
    # Create targets once for all sources
    targets = [FixedTarget(coord.SkyCoord(ra=src.ra_deg * u.deg, dec=src.dec_deg * u.deg), name=src.name) 
               for src in valid_sources]
    
    for src, target in zip(valid_sources, targets):
        if require_all_stations:
            # Check all stations
            all_visible = True
            src_min_elevs = []
            for station in stations:
                if not station.is_ever_observable(times, target):
                    all_visible = False
                    break
                elevations = station.elevation(times, target)
                if not np.all(elevations.value >= min_el_deg):
                    all_visible = False
                    break
                src_min_elevs.append(np.min(elevations.value))
            if all_visible:
                visible_sources.append(src)
                min_elevations.append(np.min(src_min_elevs))
        else:
            # Check any station
            for station in stations:
                if station.is_ever_observable(times, target):
                    elevations = station.elevation(times, target)
                    if np.any(elevations.value >= min_el_deg):
                        visible_sources.append(src)
                        min_elevations.append(np.min(elevations.value))
                        break
    
    # Sort by flux in the specified band
    visible_sources.sort(key=lambda s: s.get_flux_at_band(bands_to_check[0])[1], reverse=True)
    min_elevations = [e for _, e in sorted(zip(visible_sources, min_elevations), 
                       key=lambda x: x[0].get_flux_at_band(bands_to_check[0])[1], reverse=True)]
    
    return visible_sources, min_elevations


def get_nearby_sources(source: CalibratorSource | Source, max_separation: u.Quantity = 5.0 * u.deg, 
                       min_unresolved_flux: float = 0.0, band: Optional[str] = None, 
                       catalog: Optional[RFCCatalog] = None, 
                       n_sources: Optional[int] = None) -> list[tuple[CalibratorSource, float]]:
    """Find calibrator sources near a target source."""
    if catalog is None:
        catalog = RFCCatalog()

    if band is not None:
        if band in _WAVELENGTH_BANDS:
            band = _WAVELENGTH_BANDS[band]
        rfc_band = _RFC_BANDS.get(band.lower(), band.lower())
        bands_to_check = [rfc_band]
        ignore_min_flux = False
    else:
        bands_to_check = ['s', 'c', 'x', 'u', 'k']
        ignore_min_flux = True
    
    max_sep_deg = max_separation.to(u.deg).value if hasattr(max_separation, 'to') else max_separation
    ra_arr, dec_arr = catalog._get_coord_arrays()
    if len(ra_arr) == 0:
        return []
    
    separations_deg = _angular_separation(source.coord.ra.deg, source.coord.dec.deg, ra_arr, dec_arr)
    valid_mask = (separations_deg <= max_sep_deg) & (separations_deg > 0)
    if not ignore_min_flux:
        flux_mask = np.array([src.get_flux_at_band(bands_to_check[0])[1] > min_unresolved_flux 
                             for src in catalog.sources])
        valid_mask = valid_mask & flux_mask
    
    valid_indices = np.where(valid_mask)[0]
    nearby = [(catalog.sources[i], separations_deg[i]) for i in valid_indices]
    nearby.sort(key=lambda x: x[1])
    return nearby[:n_sources] if n_sources is not None else nearby


def main_fringe():
    usage = "%(prog)s [-h] OPTIONS"
    description = "Find fringe finder sources for VLBI observations"
    parser = argparse.ArgumentParser(description=description, prog="planobs_fringefinder", usage=usage, 
                                  formatter_class=RawTextRichHelpFormatter)
    parser.add_argument('-s', '--stations', type=str, nargs='+', required=True, 
                        help="List of antenna codenames or names that will participate in the observation.")
    parser.add_argument('-t', '--starttime', type=str, required=True, 
                        help="Start of the observation in format 'YYYY-MM-DD HH:MM' (UTC).")
    parser.add_argument('-d', '--duration', type=float, required=True, help="Duration of the observation in hours.")
    parser.add_argument('-b', '--band', type=str, required=True, 
                        help="Observing band (e.g., '18cm', '21cm', '13cm', '6cm', '5cm', '3.6cm', '2cm', '1.3cm', '0.7cm').")
    parser.add_argument('--min-flux', type=float, default=0.5, 
                        help="Minimum unresolved flux threshold in Jy (default: 0.5).")
    parser.add_argument('--min-elevation', type=float, default=20.0, 
                        help="Minimum elevation in degrees (default: 20).")
    parser.add_argument('-l', '--max-lines', type=int, default=20, 
                        help="Maximum number of sources to return (default: 20).")
    parser.add_argument('--require-all', action='store_true', default=False, 
                        help="Require source to be visible by ALL stations (default: False).")
    parser.add_argument('--station-catalog', type=str, default=None, help="Path to custom station catalog file.")
    args = parser.parse_args()
    obs._STATIONS = obs.Stations(filename=args.station_catalog)
    stations_list = []
    for s in args.stations:
        try:
            a_station = obs._STATIONS[s.strip()].codename
            if a_station not in stations_list:
                stations_list.append(a_station)
        except KeyError:
            rprint(f"[bold red]The station {s} is not known.[/bold red]")
            sys.exit(1)

    stations_obj = obs._STATIONS.filter_antennas(stations_list)
    if not stations_obj:
        rprint("[bold red]No valid antennas have been selected.[/bold red]")
        sys.exit(1)
    
    # Convert wavelength to RFC band if provided
    if args.band in _WAVELENGTH_BANDS:
        rfc_band = _WAVELENGTH_BANDS[args.band]
    else:
        rfc_band = _RFC_BANDS.get(args.band.lower(), args.band.lower())
    band_display = args.band.upper()
    band_str = f"at {band_display}"
    
    times = Time(args.starttime, scale='utc') + np.arange(0, args.duration + 0.1, 0.1) * u.hour
    sources, min_elevs = get_fringe_finder_sources(stations_obj, times, band=rfc_band, 
                                               min_elevation=args.min_elevation * u.deg, 
                                               min_flux=args.min_flux * u.Jy, 
                                               require_all_stations=args.require_all)
    if not sources:
        rprint(f"[bold red]No fringe finder candidates found {band_str} above {args.min_elevation} degrees elevation and with a unresolved flux above {args.min_flux} Jy.[/bold red]")
        sys.exit(0)

    rprint(f"\n[bold green]Found {len(sources)} fringe finder candidates {band_str} above "
           f"{args.min_elevation} degrees elevation and with a unresolved flux above {args.min_flux} Jy:[/bold green]")
    
    table = Table(show_header=True, header_style="bold", show_lines=False, box=box.SIMPLE)
    table.add_column("Name", style="", width=17)
    table.add_column("IVS Name", style="", width=10)
    table.add_column("Min elev. (deg)", justify="right", style="", width=10)
    table.add_column("Total flux (Jy)", justify="right", style="", width=12)
    table.add_column("Unresolved (Jy)", justify="right", style="", width=13)
    table.add_column("Bands", justify="right", style="", width=10)
    table.add_column("url", style="", width=10)
    
    for src, min_elev in zip(sources[:args.max_lines], min_elevs[:args.max_lines]):
        resolved_flux, unresolved_flux = src.get_flux_at_band(args.band)
        
        table.add_row(src.name, src.ivsname, f"{min_elev if min_elev > 0.0 else 0:>6.1f}",
                      f"{resolved_flux:>12.2f}", f"{unresolved_flux:>13.2f}",
                      src.get_observed_bands(), f"[link={src.get_astrogeo_link()}]AstroGeo[/link]")
    
    rprint(table)
    if len(sources) > args.max_lines:
        rprint(f"\n... and {len(sources) - args.max_lines} more sources.")
    sys.exit(0)


def main_phasecal():
    usage = "%(prog)s [-h] OPTIONS"
    description = "Find phase calibrator sources near a target source"
    parser = argparse.ArgumentParser(description=description, prog="planobs_phasecal", usage=usage, 
                                  formatter_class=RawTextRichHelpFormatter)
    parser.add_argument('-t', '--target', type=str, required=True, 
                        help="Target source name (J2000 or IVS name from RFC catalog).")
    parser.add_argument('-b', '--band', type=str, default=None, 
                        help="Observing band (e.g., '18cm', '21cm', '13cm', '6cm', '5cm', '3.6cm', '2cm', '1.3cm', '0.7cm'). If not provided, checks all bands.")
    parser.add_argument('--max-separation', type=float, default=5.0, 
                        help="Maximum angular separation in degrees (default: 5.0).")
    parser.add_argument('--min-flux', type=float, default=0.1, 
                        help="Minimum unresolved flux threshold in Jy (default: 0.1).")
    parser.add_argument('-n', '--n-sources', type=int, default=None, 
                        help="Maximum number of sources to return (default: all).")
    parser.add_argument('--catalog-file', type=str, default=None, 
                        help="Path to custom RFC catalog file.")
    args = parser.parse_args()
    
    if args.band is not None:
        if args.band not in _WAVELENGTH_BANDS and 'cm' in args.band.lower():
            original_band = args.band
            args.band = _round_to_nearest_wavelength(args.band)
            if args.band != original_band:
                rprint(f"[bold yellow]Note: Band '{original_band}' rounded to '{args.band}'[/bold yellow]")
        
        if args.band in _WAVELENGTH_BANDS:
            rfc_band = _WAVELENGTH_BANDS[args.band]
        else:
            rfc_band = _RFC_BANDS.get(args.band.lower(), args.band.lower())
        band_display = args.band.upper()
    else:
        rfc_band = None
        band_display = "all bands"
    
    catalog = RFCCatalog(catalog_filename=args.catalog_file, band='c', min_flux=0.0)
    target = catalog.get_source(args.target)
    if not target:
        try:
            target = Source.source_from_str(args.target, source_type=SourceType.TARGET)
        except Exception:
            rprint(f"[bold red]Target source '{args.target}' could not be parsed or found in the catalogs.[/bold red]")
            sys.exit(1)

    nearby = get_nearby_sources(target, max_separation=args.max_separation * u.deg, 
                               min_unresolved_flux=args.min_flux, band=rfc_band, 
                               catalog=catalog, n_sources=args.n_sources)
    if not nearby:
        band_str = f"at {band_display}"
        rprint(f"[bold red]No phase calibrator candidates found near {target.name} "
               f"({target.coord.to_string('hmsdms')}) {band_str}.[/bold red]")
        sys.exit(1)

    band_str = f"at {band_display}"
    rprint(f"\n[bold green]Found {len(nearby)} phase calibrator candidates near {target.name} "
           f"({target.coord.to_string('hmsdms')}) {band_str}:[/bold green]")
    
    table = Table(show_header=True, header_style="bold", show_lines=False, box=box.SIMPLE)
    table.add_column("Name", style="")
    table.add_column("IVS Name", style="")
    table.add_column("Separation (deg)", justify="right", style="")
    table.add_column("Total Flux (Jy)", justify="right", style="")
    table.add_column("Unresolved (Jy)", justify="right", style="")
    table.add_column("Bands", justify="right", style="", width=10)
    table.add_column("url", style="")
    
    for src, sep in nearby:
        if args.band is not None:
            resolved_flux, unresolved_flux = src.get_flux_at_band(args.band)
        else:
            resolved_flux = np.mean([src.get_flux_at_band(b)[0] for b in ['s', 'c', 'x', 'u', 'k']])
            unresolved_flux = np.mean([src.get_flux_at_band(b)[1] for b in ['s', 'c', 'x', 'u', 'k']])
        
        table.add_row(src.name, src.ivsname, f"{sep:.2f}", f"{resolved_flux:.2f}",
                      f"{unresolved_flux:.2f}", src.get_observed_bands(),
                      f"[link={src.get_astrogeo_link()}]AstroGeo[/link]")
    
    rprint(table)
    sys.exit(0)
