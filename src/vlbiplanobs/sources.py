from typing import Optional, Union, Self, Sequence
from importlib import resources
import re
import logging
import threading
import functools
import numpy as np
import tomllib
import operator
from functools import reduce
from pathlib import Path
from enum import Enum, auto
from dataclasses import dataclass
from astropy import units as u
from astropy.time import Time
from astropy import coordinates as coord
from astroplan import FixedTarget


__all__ = ['SourceNotVisible', 'Source', 'SourceType', 'Scan', 'ScanBlock']

_log = logging.getLogger(__name__)

# Module-level RFC catalog cache. Only ever assigned a fully-built dict, under _RFC_CATALOG_LOCK.
_RFC_CATALOG_CACHE: Optional[dict[str, tuple[str, str, str]]] = None
_RFC_CATALOG_LOCK = threading.Lock()

# Conservative whitelist of characters allowed in a source name sent to online resolvers (Sesame).
_VALID_SOURCE_NAME = re.compile(r"^[\w\s+\-.*:()/']{1,80}$")

"""Defines an observation, which basically consist of a given network of stations,
observing a target source for a given time range and at an observing band.
"""


class SourceNotVisible(Exception):
    """Exception produced when a given target source cannot be observed for any
    antenna in the network.
    """
    pass


@dataclass
class FluxMeasurement:
    """Stores the fluxes related to the given source at a particular band (or frequency).

    Attributes
    ----------
    resolved : u.Quantity
        Refers to the peak brightness of the source. Meaning the peak expected in a map for the source in
        a VLBI map, i.e. the unresolved flux of the source.
    unresolved : u.Quantity
        Refers to the total flux density of the source. Meaning the observed flux on the shortest baselines.
    """
    resolved: u.Quantity
    unresolved: u.Quantity


class SourceFlux(object):
    """Fluxes attributed to a given source.

    It provides the flux density (total flux) and peak flux (flux on the longest baselines)
    for a particular frequency.
    """

    def __init__(self, band_flux: dict[str, FluxMeasurement]):
        """Initializes a SourceFlux, which contains the peak flux and flux density of a given source at the given band.

        Parameters
        ----------
        band_flux : dict[str, FluxMeasurement]
            Dictionary of the form {band: FluxMeasurement}, with 'band' the different bands at which the
            flux measurements are referring to, and the corresponding FluxMeasurement object.
        """
        self._data: dict[str, FluxMeasurement] = band_flux

    def bands(self) -> tuple[str]:
        """Returns the bands at which there is flux information.

        Returns
        -------
        tuple[str]
            Tuple of strings representing the bands.
        """
        return tuple(self._data.keys())  # type: ignore

    def flux_density(self, band: str) -> FluxMeasurement:
        """Returns the flux density measurements associated to the source at the given band.

        Parameters
        ----------
        band : str
            The band at which the flux density measurements are referring to.

        Returns
        -------
        FluxMeasurement
            The flux density measurements.

        Raises
        ------
        KeyError
            If band is not available.
        """
        return self._data[band].unresolved

    def peak_flux(self, band: str) -> FluxMeasurement:
        """Returns the peak flux (brightness) measurements associated to the source at the given band.

        Parameters
        ----------
        band : str
            The band at which the flux density measurements are referring to.

        Returns
        -------
        FluxMeasurement
            The peak flux measurements.

        Raises
        ------
        KeyError
            If band is not available.
        """
        return self._data[band].resolved

    def has_band(self, band: str) -> bool:
        """Returns if the band is present in the FluxMeasurement."""
        return band in self._data

    def __contains__(self, band: str) -> bool:
        """Returns if the band is present in the FluxMeasurement."""
        return band in self._data

    def __getitem__(self, band: str) -> FluxMeasurement:
        """Returns the flux measurements associated to the source at the given band.

        Parameters
        ----------
        band : str
            The band at which the flux density measurements are referring to.

        Returns
        -------
        FluxMeasurement
            The flux density measurements.

        Raises
        ------
        KeyError
            If band is not available.
        """
        return self._data[band]

    def __setitem__(self, band: str, flux: FluxMeasurement):
        """Sets the flux measurements associated to the source at the given band.

        Parameters
        ----------
        band : str
            The band at which the flux density measurements are referring to.
        flux : FluxMeasurement
            The flux density measurements.
        """
        self._data[band] = flux

    def add_band(self, band: str, flux: FluxMeasurement):
        """Adds a new measurement of the flux at a new frequency, or overwrites a previous one if the band exists.

        Parameters
        ----------
        band : str
            The band at which the flux density measurements are referring to.
        flux : FluxMeasurement
            The flux density measurements.

        Raises
        ------
        TypeError
            If 'band' is not a str or 'flux' is not a FluxMeasurement.
        """
        if (not isinstance(band, str)) or (not isinstance(flux, FluxMeasurement)):
            raise TypeError("Expected 'band' to be a str and 'flux' to be a FluxMeasurement.")

        self._data[band] = flux


class SourceType(Enum):
    """Types of sources in a regular VLBI observation."""
    TARGET = auto()
    PHASECAL = auto()
    FRINGEFINDER = auto()
    AMPLITUDECAL = auto()
    CHECKSOURCE = auto()
    POLCAL = auto()
    PULSAR = auto()
    UNKNOWN = auto()


class Source(FixedTarget):
    """Defines a target source located at some coordinates and with a given name."""

    def __init__(self, name: str,
                 coordinates: Optional[Union[str, coord.SkyCoord]] = None,
                 source_type: SourceType = SourceType.UNKNOWN,
                 flux: Optional[SourceFlux] = None,
                 notes: Optional[str] = None,
                 other_names: Optional[list[str]] = None, **kwargs):
        """Initializes a Source object.

        Parameters
        ----------
        name : str
            Name associated to the source.
        coordinates : str or astropy.coordinates.SkyCoord, optional
            Coordinates of the target source in a str format recognized by
            astropy.coordinates.SkyCoord (e.g. XXhXXmXXs XXdXXmXXs).
            J2000 coordinates are assumed. If not provided, name must be a source name
            recognized by the RFC catalog or astroquery.
        source_type : SourceType, optional
            Defines the type of the source. Default is UNKNOWN.
        flux : SourceFlux, optional
            Estimated flux density of the source at some given frequencies.
        notes : str, optional
            Some notes that you want to add for further information on the source.
        other_names : list[str], optional
            A list of other possible names that the source may have.
        **kwargs
            Keyword arguments to be passed to astropy.coordinates.SkyCoord() if needed.
            For example, the 'unit=' parameter.

        Notes
        -----
        If both name and coordinates are provided, the given coordinates will be used for the given source.

        Raises
        ------
        NameResolveError
            If the name is not recognized (and no coordinates are provided).
        ValueError
            If the coordinates have an unrecognized format.
        AttributeError
            If neither name or coordinates are provided, or name is empty.
        """
        if not isinstance(name, str):
            raise ValueError("'name' for Source needs to be a string (single source allowed).")

        if not isinstance(source_type, SourceType):
            raise ValueError("source_type must be a SourceType value.")

        if coordinates is None:
            coordinates = self.get_coordinates_from_name(name)

        if isinstance(coordinates, coord.SkyCoord) and not kwargs:
            sky_coord = coordinates
        else:
            sky_coord = coord.SkyCoord(coordinates, **kwargs)

        super().__init__(sky_coord, name)
        self._type = source_type
        self._flux = flux
        self._notes = notes
        self._other_names = other_names if other_names is not None else list()

    # @property
    # def coordinates(self) -> coord.SkyCoord
    #     return self.
    @property
    def other_names(self) -> list[str]:
        """List of other possible names to refer to this source."""
        return self._other_names

    @other_names.setter
    def other_names(self, other_names: list[str]):
        """Sets a list of other possible names to refer to this source."""
        self._other_names = other_names

    @property
    def type(self) -> SourceType:
        """Type of the source."""
        return self._type

    @property
    def flux(self) -> Optional[SourceFlux]:
        """Estimated flux of the source."""
        return self._flux

    @property
    def notes(self) -> Optional[str]:
        """Notes on the source."""
        return self._notes

    @staticmethod
    def get_coordinates_from_name(src_name: str) -> coord.SkyCoord:
        """Returns the coordinates of a source by searching for them given the source name.

        First it searches in the RFC catalog, and if not found, in the ICRS catalogs.

        Parameters
        ----------
        src_name : str
            Name of the source to search for.

        Returns
        -------
        astropy.coordinates.SkyCoord
            The coordinates of the source, if found.

        Raises
        ------
        NameResolveError
            If there is no connection or unable to find ICRS sources.
        ValueError
            If the coordinates have an unrecognized format, or the name contains characters
            not allowed for online resolution.
        """
        try:
            return Source.get_rfc_coordinates(src_name)
        except ValueError:
            return resolve_name_online(src_name)

    @classmethod
    def source_from_name(cls, src_name: str, source_type: SourceType = SourceType.TARGET) -> Self:
        """Returns a Source object by finding the coordinates from its name.

        Parameters
        ----------
        src_name : str
            Name of the source to search for.
        source_type : SourceType, optional
            Type of the source. Default is TARGET.

        Returns
        -------
        Self
            A new Source object with coordinates from either RFC catalog or ICRS catalogs.

        Raises
        ------
        NameResolveError
            If the source cannot be found in any catalog.
        ValueError
            If the coordinates have an unrecognized format.
        """
        return cls(src_name, coordinates=Source.get_coordinates_from_name(src_name), source_type=source_type)

    @staticmethod
    def parse_source_spec(spec: str) -> tuple[Optional[str], Optional[str]]:
        """Parse a source specification that may contain 'name/coordinates'.

        Parameters
        ----------
        spec : str
            Source spec in one of these forms:
            a) 'name/coordinates' — both a custom name and explicit coordinates.
            b) 'coordinates' only (contains h/m/d/s or ':' patterns).
            c) 'name' only — to be resolved via catalog lookup.

        Returns
        -------
        tuple[str | None, str | None]
            (name, coord_str). At least one will be non-None.
            If '/' is present, name is the part before '/' and coord_str the part after.
            If no '/', returns (None, spec) when spec looks like coordinates,
            or (spec, None) when spec looks like a name.
        """
        if '/' in spec:
            name, coord_str = spec.split('/', 1)
            return name.strip(), coord_str.strip()

        # No '/' — check if it looks like coordinates
        if all(char in spec for char in ('h', 'm', 'd', 's')) or ':' in spec:
            return None, spec

        return spec, None

    @classmethod
    def _parse_coord_str(cls, coord_str: str) -> coord.SkyCoord:
        """Parse a coordinate string in 'XXhXXmXXs XXdXXmXXs' or 'HH:MM:SS DD:MM:SS' format.

        Parameters
        ----------
        coord_str : str
            Coordinate string to parse.

        Returns
        -------
        astropy.coordinates.SkyCoord

        Raises
        ------
        ValueError
            If the coordinate format is invalid.
        """
        if all(char in coord_str for char in ('h', 'm', 'd', 's')):
            return coord.SkyCoord(coord_str)

        if ':' in coord_str:
            temp = coord_str
            for char in ('h', 'm', 'd', 'm'):
                temp = temp.replace(':', char, 1)
            return coord.SkyCoord(temp)

        raise ValueError(f"Cannot parse coordinate string: '{coord_str}'")

    @classmethod
    def source_from_str(cls, src: str, source_type: SourceType = SourceType.TARGET) -> Self:
        """Returns a Source object from a source name, coordinate string, or 'name/coordinates'.

        Parameters
        ----------
        src : str
            One of:
            a) 'name/coordinates' — use the given name with explicit coordinates.
            b) Coordinates in 'XXhXXmXXs XXdXXmXXs' or 'HH:MM:SS DD:MM:SS' format.
            c) A source name to be looked up in catalogs.
        source_type : SourceType, optional
            Type of the source. Default is TARGET.

        Returns
        -------
        Self
            A new Source object with the specified coordinates.

        Raises
        ------
        ValueError
            If the coordinate format is invalid.
        NameResolveError
            If the source name cannot be found in catalogs.
        """
        name, coord_str = cls.parse_source_spec(src)

        if name is not None and coord_str is not None:
            return cls(name, coordinates=cls._parse_coord_str(coord_str), source_type=source_type)

        if coord_str is not None:
            return cls('target', coordinates=cls._parse_coord_str(coord_str), source_type=source_type)

        return cls.source_from_name(name, source_type)

    @staticmethod
    def get_rfc_coordinates(src_name: str) -> coord.SkyCoord:
        """Returns the coordinates of the object by searching the provided name through the RFC catalog.

        This method uses an in-memory cache of the RFC catalog for fast lookups,
        avoiding subprocess calls and file I/O on repeated accesses.

        Parameters
        ----------
        src_name : str
            Name of the source to search for in the RFC catalog (case-insensitive).

        Returns
        -------
        astropy.coordinates.SkyCoord
            The coordinates of the source from the RFC catalog.

        Raises
        ------
        ValueError
            If the source name is not found in the RFC catalog.
        RuntimeError
            If no RFC catalog files are found.
        """
        catalog = _load_rfc_catalog()
        src_upper = src_name.upper()

        if src_upper in catalog:
            _, _, coord_str = catalog[src_upper]
            return coord.SkyCoord(coord_str)

        raise ValueError(f"The source {src_name} was not found in the RFC catalog.")

    def sun_separation(self, times: Time) -> Sequence[u.Quantity]:
        """Returns the separation of the source to the Sun at the given epoch(s).

        Parameters
        ----------
        times : astropy.time.Time
            An array of times defining the duration of the observation. The first time
            defines the start and the last one the end of the observation. Higher time
            resolution provides more precise values but increases computation time.

        Returns
        -------
        Sequence[astropy.units.Quantity]
            The angular separation between the source and the Sun at each given time.
        """
        return self.coord.transform_to(coord.GCRS(obstime=times)).separation(coord.get_sun(times))

    def sun_constraint(self, min_separation: u.Quantity, times: Optional[Time] = None) -> Time:
        """Returns times when the Sun is too close to observe the source.

        Parameters
        ----------
        min_separation : astropy.units.Quantity
            Minimum allowed angular separation between source and Sun.
            See `freqsetups.solar_separations` for default values per band.
        times : astropy.time.Time, optional
            Times to check Sun separation. If None, checks full current year
            at 1-day resolution.

        Returns
        -------
        astropy.time.Time
            Times when Sun is closer than min_separation. Empty if never too close.
        """
        if times is None:
            times = Time(f"{Time.now().datetime.year}-01-01") + np.arange(0, 365, 1)*u.day

        sun_separation = self.sun_separation(times=times)
        # if isinstance(sun_separation, list):
        #     return [times[sun_separation[i] < min_separation] for i in range(len(self.coord))]
        # else:
        return times[sun_separation < min_separation]


def _load_rfc_catalog() -> dict[str, tuple[str, str, str]]:
    """Load the RFC catalog file into memory and cache it.

    Returns
    -------
    dict[str, tuple[str, str, str]]
        Dictionary mapping source names (uppercase) to tuples of (IVS name, J2000 name, coordinate string).
        The coordinate string is in format 'XXhXXmXXs XXdXXmXXs'.

    Raises
    ------
    RuntimeError
        If no RFC catalog files are found.
    """
    global _RFC_CATALOG_CACHE

    if _RFC_CATALOG_CACHE is not None:
        return _RFC_CATALOG_CACHE

    with _RFC_CATALOG_LOCK:
        # Re-check: another thread may have finished loading while we waited for the lock.
        if _RFC_CATALOG_CACHE is not None:
            return _RFC_CATALOG_CACHE

        rfc_files = tuple(r.name for r in resources.files("vlbiplanobs.data").iterdir()
                          if r.is_file() and 'rfc' in r.name)
        if not rfc_files:
            raise RuntimeError("No RFC files found under the 'data' folder.")

        # Build locally; publish only once complete so no partial cache is ever visible.
        catalog: dict[str, tuple[str, str, str]] = {}
        with resources.as_file(resources.files("vlbiplanobs.data").joinpath(sorted(rfc_files)[-1])) as rfcfile:
            with open(rfcfile, 'r') as f:
                for line in f:
                    parts = line.split()
                    if len(parts) >= 9:  # IVS name J2000name h m s d m s
                        ivs_name = parts[0]
                        j2000_name = parts[2]
                        coord_str = f"{parts[3]}h{parts[4]}m{parts[5]}s {parts[6]}d{parts[7]}m{parts[8]}s"
                        catalog[j2000_name.upper()] = (ivs_name, j2000_name, coord_str)
                        catalog[ivs_name.upper()] = (ivs_name, j2000_name, coord_str)

        _RFC_CATALOG_CACHE = catalog
        _log.info("Loaded RFC name catalog (%d keys)", len(catalog))

    return _RFC_CATALOG_CACHE


def validate_source_name(src_name: str) -> str:
    """Checks that a source name is safe to send to an online name resolver.

    Parameters
    ----------
    src_name : str
        Source name to validate.

    Returns
    -------
    str
        The stripped source name.

    Raises
    ------
    ValueError
        If the name is not a str, is empty, longer than 80 characters, or contains characters
        other than letters, digits, whitespace and + - . * : ( ) / '.
    """
    if not isinstance(src_name, str):
        raise ValueError(f"Source name must be a str, got {type(src_name).__name__}.")

    name = src_name.strip()
    if not _VALID_SOURCE_NAME.match(name):
        raise ValueError(f"Invalid source name {src_name!r}: it must be 1-80 characters long and contain only "
                         "letters, digits, spaces and + - . * : ( ) / '.")

    return name


@functools.lru_cache(maxsize=1024)
def _cached_icrs_coordinates(name: str) -> coord.SkyCoord:
    """Online (Sesame) name resolution, cached per name. Failures raise and are therefore not cached."""
    _log.info("Resolving source name '%s' via online services", name)
    return coord.get_icrs_coordinates(name)


def resolve_name_online(src_name: str) -> coord.SkyCoord:
    """Resolves a source name into coordinates using online services (Sesame: SIMBAD/NED/VizieR).

    Names are validated before any network request, and successful results are cached (up to 1024 names).

    Parameters
    ----------
    src_name : str
        Source name to resolve.

    Returns
    -------
    astropy.coordinates.SkyCoord
        ICRS coordinates of the source.

    Raises
    ------
    ValueError
        If the name fails validation (see `validate_source_name`).
    astropy.coordinates.name_resolve.NameResolveError
        If the name cannot be resolved (or no connection is available).
    """
    return _cached_icrs_coordinates(validate_source_name(src_name))


def _format_dec_coordinate(dec_str: str) -> str:
    """Ensure declination string has proper sign prefix for SkyCoord parsing.

    Parameters
    ----------
    dec_str : str
        Declination string in format 'DD:MM:SS' or '+DD:MM:SS' or '-DD:MM:SS'.

    Returns
    -------
    str
        Formatted declination string with guaranteed sign prefix.
    """
    if not dec_str.startswith(('+', '-')):
        dec_str = '+' + dec_str
    return dec_str


def _resolve_coord_str(entry: dict, name: str) -> Optional[str]:
    """Return a SkyCoord-parseable coordinate string for a catalog entry.

    Tries, in order:
    1. 'coordinates' key in the TOML entry (explicit RA/Dec).
    2. RFC catalog lookup by name.
    3. Online name resolution (SIMBAD/NED/VizieR).

    Returns None and logs a warning if all methods fail.

    Parameters
    ----------
    entry : dict
        Parsed TOML sub-dict for a phasecal, checksource, or target.
    name : str
        Source name, used for catalog/online lookup when coordinates are absent.

    Returns
    -------
    str or None
        SkyCoord-parseable coordinate string, or None if resolution fails.
    """
    if 'coordinates' in entry:
        c = entry['coordinates']
        ra = c['RA'].replace(':', 'h', 1).replace(':', 'm') + 's'
        dec = _format_dec_coordinate(c['Dec']).replace(':', 'd', 1).replace(':', 'm') + 's'
        return f"{ra} {dec}"
    # Try RFC catalog
    try:
        sky = Source.get_rfc_coordinates(name)
        return sky.to_string('hmsdms')
    except (ValueError, RuntimeError):
        pass
    # Try online resolution
    try:
        sky = resolve_name_online(name)
        return sky.to_string('hmsdms')
    except (ValueError, OSError, coord.name_resolve.NameResolveError) as e:
        _log.warning("Online resolution failed for '%s': %s", name, e)
    _log.warning(
        "Could not resolve coordinates for '%s': not in catalog, RFC, or online resolvers. "
        "Skipping this source.", name)
    return None


def _scan_from_toml(entry: dict, source: 'Source', default_every: int) -> 'Scan':
    """Creates a Scan from a TOML catalog entry.

    If the entry has no 'duration', the Scan default duration is used (the kwarg is omitted, never None).

    Parameters
    ----------
    entry : dict
        Parsed TOML sub-dict (target, phasecal or checksource) with optional 'duration' (min) and 'every'.
    source : Source
        The source observed by the scan.
    default_every : int
        Value of 'every' to use when the entry does not define it.

    Returns
    -------
    Scan
    """
    kwargs = {'every': int(entry['every']) if 'every' in entry else default_every}
    if 'duration' in entry:
        kwargs['duration'] = float(entry['duration'])*u.min

    return Scan(source=source, **kwargs)


def _scans_from_toml_entry(src: dict, name: str, main_type: 'SourceType') -> list['Scan']:
    """Builds the scans (phasecal, checksource, main source) for one [[target]]/[[pulsar]] TOML entry.

    Parameters
    ----------
    src : dict
        Parsed TOML entry for the main source, with optional 'phasecal' and 'checksource' sub-tables.
    name : str
        Name of the main source.
    main_type : SourceType
        Type of the main source (TARGET or PULSAR).

    Returns
    -------
    list[Scan]
        Scans in order: phasecal (if resolvable), checksource (if resolvable), main source.

    Raises
    ------
    ValueError
        If the coordinates of the main source cannot be resolved.
    """
    scans = []
    for sub_key, sub_type, default_every in (('phasecal', SourceType.PHASECAL, -1),
                                             ('checksource', SourceType.CHECKSOURCE, 4)):
        sub = src.get(sub_key)
        if not sub or not sub.get('name'):
            continue

        sub_coord_str = _resolve_coord_str(sub, sub['name'])
        if sub_coord_str is not None:
            sub_source = Source(name=sub['name'], coordinates=sub_coord_str, source_type=sub_type)
            scans.append(_scan_from_toml(sub, sub_source, default_every))

    # Main source (coordinates required for the primary target)
    src_coord_str = _resolve_coord_str(src, name)
    if src_coord_str is None:
        raise ValueError(f"Could not resolve coordinates for {main_type.name.lower()} '{name}'. "
                         "Add 'coordinates' to its catalog entry.")

    main_source = Source(name=name, coordinates=src_coord_str, source_type=main_type)
    scans.append(_scan_from_toml(src, main_source, -1))
    return scans


class SourceCatalog:
    """Catalog of source blocks for scheduling."""

    def __init__(self, personal_catalog: Optional[str] = None):
        """Initializes a SourceCatalog.

        Parameters
        ----------
        personal_catalog : str, optional
            Path to a personal TOML catalog file to read.
        """
        self._blocks: dict[str, dict[str, ScanBlock]] = dict()
        # Per-instance caches for source_names()/sources(), keyed by include_calibrators.
        # Invalidated whenever a catalog file is (re)read.
        self._cache_source_names: dict[bool, list[str]] = {}
        self._cache_sources: dict[bool, dict[str, Source]] = {}
        if personal_catalog is not None:
            self.read_personal_catalog(personal_catalog)

    @property
    def blocknames(self):
        """List of all block names in the catalog."""
        return [bb for b in self._blocks.values() for bb in b.keys()]

    @property
    def blocks(self):
        """Dictionary of all blocks in the catalog."""
        return {bb_key: bb_value for b in self._blocks.values() for bb_key, bb_value in b.items()}

    @property
    def targets(self):
        """Target blocks in the catalog."""
        return self._blocks.get('targets', {})

    @property
    def pulsars(self):
        """Pulsar blocks in the catalog, or None if none."""
        return self._blocks['pulsars'] if 'pulsars' in self._blocks else None

    @property
    def ampcals(self):
        """Amplitude calibrator blocks in the catalog, or None if none."""
        return self._blocks['ampcals'] if 'ampcals' in self._blocks else None

    @property
    def fringefinders(self):
        """Fringe finder blocks in the catalog, or None if none."""
        return self._blocks['fringefinders'] if 'fringefinders' in self._blocks else None

    @property
    def polcals(self):
        """Polarization calibrator blocks in the catalog, or None if none."""
        return self._blocks['polcals'] if 'polcals' in self._blocks else None

    def source_names(self, include_calibrators: bool = False) -> list[str]:
        """Returns the names of all sources in the database.

        Parameters
        ----------
        include_calibrators : bool, optional
            If True, includes calibrator sources. Default is False.

        Returns
        -------
        list[str]
            List of source names.
        """
        if include_calibrators not in self._cache_source_names:
            if include_calibrators:
                self._cache_source_names[include_calibrators] = \
                    [s.name for b in self._blocks.values() for bs in b.values() for s in bs.sources()]
            else:
                self._cache_source_names[include_calibrators] = \
                    [s.name for b in self._blocks.get('targets', {}).values() for s in b.sources()]

        return self._cache_source_names[include_calibrators]

    def sources(self, include_calibrators: bool = False) -> dict[str, Source]:
        """Returns all sources.

        Parameters
        ----------
        include_calibrators : bool, optional
            If True, includes calibrator sources. Default is False.

        Returns
        -------
        dict[str, Source]
            Dictionary mapping source names to Source objects.
        """
        if include_calibrators not in self._cache_sources:
            if include_calibrators:
                self._cache_sources[include_calibrators] = \
                    {s.name: s for b in self._blocks.values() for bs in b.values() for s in bs.sources()}
            else:
                self._cache_sources[include_calibrators] = \
                    {s.name: s for b in self._blocks.get('targets', {}).values() for s in b.sources()}

        return self._cache_sources[include_calibrators]

    def __contains__(self, item: str):
        """Returns True if the item is in the catalog."""
        return item in self.blocknames

    def __getitem__(self, item: str):
        """Returns the block with the given name."""
        for key in self._blocks:
            if item in self._blocks[key]:
                return self._blocks[key][item]

        raise KeyError(f"The item {item} is not in the catalog.")

    def read_personal_catalog(self, path: str):
        """Reads a TOML file containing source information for scheduling.

        Parameters
        ----------
        path : str
            Path to the TOML containing the source catalog.

        Raises
        ------
        FileNotFoundError
            If the catalog file cannot be found.
        tomllib.TOMLDecodeError
            If the TOML file is malformed.
        ValueError
            If the source type is not recognized.
        """
        self._cache_source_names.clear()
        self._cache_sources.clear()
        with open(path, 'rb') as sources_toml:
            catalog = tomllib.load(sources_toml)

            for toml_key, block_key, main_type in (('pulsar', 'pulsars', SourceType.PULSAR),
                                                   ('target', 'targets', SourceType.TARGET)):
                if toml_key not in catalog:
                    continue

                if block_key not in self._blocks:
                    self._blocks[block_key] = dict()

                for src in catalog[toml_key]:
                    name = src.get('name', 'unknown')
                    self._blocks[block_key][name] = ScanBlock(_scans_from_toml_entry(src, name, main_type))

    def read_rfc_catalog(self, path: Optional[Union[str, Path]] = None):
        """Reads the RFC catalog file.

        Parameters
        ----------
        path : str or Path, optional
            Path to the RFC catalog file. If None, uses the default catalog.

        Raises
        ------
        FileNotFoundError
            If the catalog file cannot be found.
        ValueError
            If the catalog format is invalid.

        Notes
        -----
        This method is not implemented.
        """
        raise NotImplementedError
        # TODO: convert this to another module and use duckDB, should be much faster
        if path is None:
            path = resources.as_file(resources.files("vlbiplanobs.data").joinpath("rfc_2021_cat.txt"))

        with open(path, 'rt') as stations_catalog_path:
            pass

        # This is old code

        if path is None:
            with resources.as_file(resources.files("vlbiplanobs.data").joinpath("rfc_2021_cat.txt")) \
                                                                              as stations_catalog_path:
                all_lines = open(stations_catalog_path, 'rt').readlines()
            # stations_catalog_path = 'data/rfc_2021c_cat.txt'
        else:
            with open(path, 'rt') as stations_catalog_path:
                all_lines = stations_catalog_path.readlines()

        for aline in all_lines:
            if aline.strip()[0] == '#':
                continue

            pars = aline.strip().split()
            if len(pars) != 25:
                raise ValueError(f"Expected 25 elements in a row but {len(pars)} found: {pars}")

            a_flux = {}
            for cm, apar in zip(('13', '6', '3.6', '2', '1.3'), (14, 16, 18, 20, 22)):
                if pars[apar] != '-1.00' and pars[apar].isnumeric():
                    a_flux[cm] = float(pars[apar])*u.Jy

            # if len(a_flux) > 0:
            #     grade = 9 if all([af > 0.7*u.Jy for af in a_flux.values()]) else 5
            # else:
            #     grade = 3

            self.add(Source(name=pars[2], grade=8 if pars[0] == 'C' else 5 if pars[0] == 'N' else 3,
                            other_names=[pars[1]],
                            coordinates="{0}h{1}m{2}s {3}d{4}m{5}s".format(*pars[3:9]),
                            source_type=SourceType.PHASECAL, flux=a_flux, notes="RFC Source"),
                     label='catalog')


@dataclass
class Scan:
    """Defines a single scan of a source.

    Attributes
    ----------
    source : Source
        The source to observe.
    duration : u.Quantity
        Duration of the scan. Default is 10 minutes.
    every : int
        If positive, repeat this scan every N cycles. If -1, observe on every cycle. Default is -1.

    Raises
    ------
    ValueError
        If 'every' is 0 or lower than -1.
    """
    source: Source
    duration: u.Quantity = 10*u.min
    every: int = -1

    def __post_init__(self):
        """Validates the scan parameters."""
        if self.every == 0 or self.every < -1:
            raise ValueError(f"Scan of '{self.source.name}': 'every' must be -1 (every cycle) or a positive "
                             f"integer (every N cycles), got {self.every}.")


class ScanBlock:
    """Defines a list of scans, each of them defined as a pointing to a given source during a given time.

    A block can consist of a single target scan (for non phase-referencing observations),
    which will be repeated until the maximum observing time is filled.
    In phase-referencing observations, a scan block can consist of one or multiple target scans,
    the associated phase-referencing scans, and possible check sources to be observed every
    certain number of target scans.
    """

    def __init__(self, scans: list[Scan]):
        """Creates a block of scans.

        Ideally a block to be observed with the target scans, phase-reference calibrator source
        (if needed), and check sources.

        Parameters
        ----------
        scans : list[Scan]
            List of scans to include in the block.

        Raises
        ------
        ValueError
            If the scan block contains an empty list of scans or if any element is not a Scan.
        """
        if not scans:
            raise ValueError("The scan block cannot contain an empty list of scans.")

        if not all([isinstance(s, Scan) for s in scans]):
            raise ValueError("All elements in the list of scans must be of Scan type.")

        self._scans = scans

    @property
    def scans(self) -> list[Scan]:
        """List of scans in the block."""
        return self._scans

    def has(self, source_type: SourceType) -> bool:
        """Returns if the given source type is included among the ones observed in the provided list of scans.

        Parameters
        ----------
        source_type : SourceType
            The source type to check for.

        Returns
        -------
        bool
            True if the source type is present in the block.
        """
        return any([s.source.type is source_type for s in self.scans])

    def sources(self, source_type: Optional[SourceType] = None) -> list[Source]:
        """Returns the sources with the given source types in this block.

        Parameters
        ----------
        source_type : SourceType, optional
            The source type to filter by. If None, returns all sources.

        Returns
        -------
        list[Source]
            List of sources matching the given type, or all sources if type is None.
        """
        if source_type is None:
            return [s.source for s in self.scans]

        return [s.source for s in self.scans if s.source.type is source_type]

    def sourcenames(self, source_type: Optional[SourceType] = None) -> list[str]:
        """Returns the source names with the given source types in this block.

        Parameters
        ----------
        source_type : SourceType, optional
            The source type to filter by. If None, returns all source names.

        Returns
        -------
        list[str]
            List of source names matching the given type, or all names if type is None.
        """
        if source_type is None:
            return [s.source.name for s in self.scans]

        return [s.source.name for s in self.scans if s.source.type is source_type]

    def scan_with_sourcename(self, source_name: str) -> Optional[Scan]:
        """Returns the scan with the given source name.

        Parameters
        ----------
        source_name : str
            The source name to search for.

        Returns
        -------
        Scan or None
            The scan with the given source name, or None if not found.

        Raises
        ------
        ValueError
            If the source is not present in any scan.
        """
        for scan in self._scans:
            if scan.source.name == source_name:
                return scan

        raise ValueError(f"The source {source_name} is not present in any scan.")

    def scans_with_sources(self, source_type: SourceType) -> list[Scan]:
        """Returns the scans with the given source types in this block.

        Parameters
        ----------
        source_type : SourceType
            The source type to filter by.

        Returns
        -------
        list[Scan]
            List of scans matching the given source type.
        """
        return [s for s in self._scans if s.source.type == source_type]

    def fractional_time(self) -> dict[str, float]:
        """Returns the fractional time dedicated to each source observed within the scan block.

        Returns
        -------
        dict[str, float]
            The fraction of time estimated to be spent on each particular source assuming
            the durations of the other scans to be observed in this block. The keys of the dict
            are the source names, and the values are the fraction of time, from the total scan block
            time, spent on the source.
        """
        # Computed on every call (cheap) so that changes to the scans are always reflected.
        frac_time: dict[str, float] = {}
        # Get scans with valid durations
        valid_scans = [s for s in self.scans if s.duration is not None]

        if not valid_scans:
            return frac_time

        total_duration = sum([s.duration for s in valid_scans if s.every <= 0])

        # Over mcm_every cycles, a scan with 'every = N > 0' is observed mcm_every/N times,
        # while the base scans (every <= 0) are observed on every cycle.
        positive_every = [s.every for s in valid_scans if s.every > 0]
        if positive_every:
            mcm_every = np.lcm.reduce(positive_every)
            total_duration = total_duration*mcm_every + \
                sum([s.duration*(mcm_every/s.every) for s in valid_scans if s.every > 0])
        else:
            mcm_every = 1

        for ascan in valid_scans:
            if ascan.duration is not None:
                frac_time[ascan.source.name] = ascan.duration*mcm_every / \
                                               (ascan.every if ascan.every > 0 else 1) / total_duration

        return frac_time

    def _scans_in_cycle(self, n_loop: int) -> list[Scan]:
        """Returns the non-phasecal scans to observe in the given (1-based) cycle.

        Scans with 'every = N > 0' are observed on cycles multiple of N, replacing the target scans of that
        cycle. On any other cycle, the target scans are observed.

        Parameters
        ----------
        n_loop : int
            Cycle number, starting at 1.

        Returns
        -------
        list[Scan]
            Scans to observe in this cycle (excluding the bracketing phasecal scans).
        """
        periodic = [s for s in self.scans if s.every > 0 and n_loop % s.every == 0]
        return periodic if periodic else self.scans_with_sources(SourceType.TARGET)

    def fill(self, max_duration: u.Quantity) -> list[Scan]:
        """Given the list of scans, returns the final arrangement of scans that fills the available time.

        This will follow the following conditions:
        - Repeats the target scan within the given time. If a `phasecal` is provided, then it will always
          bracket each target scan with this `phasecal`. If two or more phasecal are provided (e.g. P1, P2),
          then it will assume a multi phase-referencing technique for the target T:
            P1 P2 T P1 P2 T P1 P2...
        - If some sources have the 'every = N > 0' constraint, then these will be scheduled every N cycles.
          For example, if no phasecal are provided and N = 3 for the C source, then: T T C T T C ...
          And if one phasecal is provided: P T P T P C P T P...

        Parameters
        ----------
        max_duration : u.Quantity
            Maximum duration to fill with scans.

        Returns
        -------
        list[Scan]
            The final arrangement of scans that fills the available time.

        Raises
        ------
        ValueError
            If the max_duration is shorter than the time of all single scans, if phase calibrator
            scans are provided without target scans, if there are no target scans, or if the
            phasecal + target cycle has zero duration.

        Notes
        -----
        Block scans should be easy! This program is not prepared for the situation when you have
        multiple targets with multiple phase reference sources mixed.
        """
        # safety Checks
        if reduce(operator.add, [s.duration.to(u.min) for s in self.scans]) > max_duration:
            raise ValueError("The max_duration of the block cannot be shorter than the time of "
                             "all single scans.")

        phasecals = self.scans_with_sources(SourceType.PHASECAL)
        targets = self.scans_with_sources(SourceType.TARGET)
        if phasecals and not targets:
            raise ValueError("If phase calibrator scans provided, then target scans must also be provided.")

        if not targets:
            raise ValueError("The scan block needs at least one target scan to be filled.")

        loop_duration = sum([s.duration.to(u.min) for s in phasecals + targets], 0*u.min)
        if loop_duration <= 0*u.min:
            raise ValueError("The phasecal + target scans of the block have zero total duration; "
                             "cannot fill the block.")

        # The phase-referencing loop needs to be closed at the end with the phasecal scans.
        last_duration = sum([s.duration.to(u.min) for s in phasecals], 0*u.min)
        main_loop: list[Scan] = []
        has_periodic = any([s.every > 0 for s in self.scans
                            if s.source.type not in (SourceType.TARGET, SourceType.PHASECAL)])
        if has_periodic:
            # As there can be multiple sources to be observed every certain scans, do it incrementally.
            # Cycles repeat with period lcm(every); a zero-duration period would loop forever.
            period = int(np.lcm.reduce([s.every for s in self.scans if s.every > 0]))
            period_duration = sum([a.duration.to(u.min) for n in range(1, period + 1)
                                   for a in phasecals + self._scans_in_cycle(n)], 0*u.min)
            if period_duration <= 0*u.min:
                raise ValueError("The scans of the block have zero total duration over a full cycle; "
                                 "cannot fill the block.")

            booked_time, n_loop = 0*u.min, 1
            while True:
                to_append = phasecals + self._scans_in_cycle(n_loop)
                cycle_duration = sum([a.duration.to(u.min) for a in to_append], 0*u.min)
                n_loop += 1
                if cycle_duration + booked_time > max_duration - last_duration:
                    break

                main_loop += to_append
                booked_time += cycle_duration
        else:
            # n = floor((max_duration - last_duration) / loop_duration)
            # ensures n * loop_duration + last_duration <= max_duration
            n_reps = int((max_duration - last_duration).to(u.min).value // loop_duration.to(u.min).value)
            main_loop += (phasecals + targets) * n_reps

        main_loop += phasecals
        return main_loop

    def __iter__(self):
        """Iterate over scans in the block."""
        yield from self._scans

    def __contains__(self, a_source_name: str):
        """Returns True if the source name is in the block."""
        return a_source_name in self.sourcenames()
