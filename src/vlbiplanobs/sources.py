from typing import Optional, Union, Self, Sequence
from importlib import resources
import subprocess
import functools
# from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import numpy as np
import yaml                # type: ignore
# from rich import print as rprint
import operator
from functools import reduce
import datetime as dt
from pathlib import Path
from enum import Enum, auto
from dataclasses import dataclass
from astropy import units as u
from astropy.time import Time
from astropy import coordinates as coord
# import polars as pl
from astroplan import FixedTarget


__all__ = ['SourceNotVisible', 'Source', 'SourceType', 'Scan', 'ScanBlock']

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

    - peak_flux : u.Quantity
        Refers to the peak brightness of the source. Meaning the peak expected in a map for the source in
        a VLBI map, i.e. the unresolved flux of the source.

    - flux_density : u.Quantity
        Refers to the total flux density of the source. Meaning the observed flux on the shortest baselines.
    """
    peak_flux: u.Quantity
    flux_density: u.Quantity


class SourceFlux(object):
    """Fluxes attributed to a given source.
    It provides the flux density (total flux) and peak flux (flux on the longest baselines)
    for a particular frequency.
    """

    def __init__(self, band_flux: dict[str, FluxMeasurement]):
        """Initializes a SourceFlux, which contains the peak flux and flux density of a given source at the
        given band.

        Input
            - band_flux : dict[str, FluxDensity]
                Dictionary of the form {band: FluxDensity}, with 'band' the different bands at which the
                flux measurements are referring to, and the corresponding FluxDensity object.
        """
        self._data: dict[str, FluxMeasurement] = band_flux

    def bands(self) -> tuple[str]:
        """Returns the bands at which there is flux information.
        """
        return tuple(self._data.keys())  # type: ignore

    def flux_density(self, band: str) -> FluxMeasurement:
        """Returns the flux density measurements associated to the source at the given band.

        Input
        - band : str
            The band at which the flux density measurements are referring to.

        Raises
            KeyError : if band is not available.
        """
        return self._data[band].flux_density

    def peak_flux(self, band: str) -> FluxMeasurement:
        """Returns the peak flux (brightness) measurements associated to the source at the given band.

        Input
        - band : str
            The band at which the flux density measurements are referring to.

        Raises
            KeyError : if band is not available.
        """
        return self._data[band].peak_flux

    def has_band(self, band: str) -> bool:
        """Returns if the band is present in the FluxMeasurement.
        """
        return band in self._data

    def __containts__(self, band: str) -> bool:
        return band in self._data

    def __getitem__(self, band: str) -> FluxMeasurement:
        return self._data[band]

    def __setitem__(self, band: str, flux: FluxMeasurement):
        self._data[band] = flux

    def add_band(self, band: str, flux: FluxMeasurement):
        """Adds a new measurement of the flux at a new frequency, or overwrites a previous
        one if the band exists.
        """
        if (not isinstance(band, str)) or (not isinstance(flux, FluxMeasurement)):
            raise TypeError("Expected 'band' to be a str and 'flux' to be a FluxMeasurement.")

        self._data[band] = flux


class SourceType(Enum):
    """Types of sources in a regular VLBI observation
    """
    TARGET = auto()
    PHASECAL = auto()
    FRINGEFINDER = auto()
    AMPLITUDECAL = auto()
    CHECKSOURCE = auto()
    POLCAL = auto()
    PULSAR = auto()
    UNKNOWN = auto()


class Source(FixedTarget):
    """Defines a target source located at some coordinates and with a given name.
    """

    def __init__(self, name: str,
                 coordinates: Optional[Union[str, coord.SkyCoord]] = None,
                 source_type: SourceType = SourceType.UNKNOWN,
                 flux: Optional[SourceFlux] = None,
                 notes: Optional[str] = None,
                 other_names: Optional[list[str]] = None, **kwargs):
        """Initializes a Source object.

        Inputs
        - name : str
            Name associated to the source.
        - coordinates : str or astropy.coordinates.SkyCoord [OPTIONAL]
            Coordinates of the target source in a str format recognized by
            astropy.coordinates.SkyCoord (e.g. XXhXXmXXs XXdXXmXXs).
            J2000 coordinates are assumed.
            If not provided, name must be a source name recognized by the RFC catalog or astroquery.
        - source_type : SourceType  [default = UNKNOWN]
            Defines the type of the source.
        - flux : SourceFlux  [default = None]
            Estimated flux density of the source at some given frequencies.
        - notes : str  [default = None]
            Some notes that you want to add for further information on the source.
        - other_names : list[str]  [default = None]
            A list of other possible names that the source may have.
        - kwargs
            keyword arguments to be passed to astropy.coordinates.SkyCoord() if needed.
            For example, the 'unit=' parameter.

        If both provided, the given coordinates will be used for the given source.

        It may raise:
        - NameResolveError: if the name is not recognized (and no coordinates are provided)
        - ValueError: if the coordinates have an unrecognized format.
        - AttributeError: if neither name or coordinates are provided, or name is empty.
        """
        if not isinstance(name, str):
            raise ValueError("'name' for Source needs to be a string (single source allowed).")

        if not isinstance(source_type, SourceType):
            raise ValueError("source_type must be a SourceType value.")

        if coordinates is None:
            coordinates = self.get_coordinates_from_name(name)

        super().__init__(coord.SkyCoord(coordinates, **kwargs), name)
        self._type = source_type
        self._flux = flux
        self._notes = notes
        self._other_names = other_names if other_names is not None else list()

    # @property
    # def coordinates(self) -> coord.SkyCoord
    #     return self.

    @property
    def other_names(self) -> list[str]:
        """Returns a list of other possible names to refer to this source.
        It may be None, if there are no other names.
        """
        return self._other_names

    @other_names.setter
    def other_names(self, other_names: list[str]):
        self._other_names = other_names

    @property
    def type(self) -> SourceType:
        return self._type

    @property
    def flux(self) -> Optional[SourceFlux]:
        """Estimated flux of the source.
        """
        return self._flux

    @property
    def notes(self) -> Optional[str]:
        """Some notes on the source
        """
        return self._notes

    @staticmethod
    def get_coordinates_from_name(src_name: str) -> coord.SkyCoord:
        """Returns the coordinates of a source by searching for them given the source name.
        First it searches in the RFC catalog, and if not found, in the ICRS catalogs.

        Inputs
        - src_name : str
            Name of the source to search for.

        Returns
        - src_coord : astropy.coordinates.SkyCoord
            The coordinates of the source, if found.

        It may raise
        - NameResolveError - if there is no connection or unable to find ICRS sources.
        - ValueError: if the coordinates have an unrecognized format.
        """
        try:
            return Source.get_rfc_coordinates(src_name)
        except ValueError:
            return coord.get_icrs_coordinates(src_name)

    @classmethod
    def source_from_name(cls, src_name: str) -> Self:
        """Returns a Source object by finding the coordinates from its name.
        It first searches wihtin the RFC catalog, and if not found, through the
        'astropy.coordinates.get_icrs_coordinates' function through the RFC catalogs.
        """
        return cls(src_name, coordinates=Source.get_coordinates_from_name(src_name))

    @classmethod
    def source_from_str(cls, src: str) -> Self:
        """Returns a Source object from the given string. This can be either the source name, and then it
        will retrieve the coordinates from the SIMBAD/NED/VizieR databases, if it's known, or the source
        coordinates (in either 'XXhXXmXXs XXdXXmXXs' or 'HH:MM:SS DD:MM:SS' formats).
        """
        if all([char in src for char in ('h', 'm', 'd', 's')]):
            return cls('target', coordinates=coord.SkyCoord(src), source_type=SourceType.TARGET)
        elif ':' in src:
            for char in ('h', 'm', 'd', 'm'):
                src = src.replace(':', char, 1)

            return cls('target', coordinates=coord.SkyCoord(src), source_type=SourceType.TARGET)
        else:
            return cls.source_from_name(src)

    @staticmethod
    def get_rfc_coordinates(src_name: str) -> coord.SkyCoord:
        """Returns the coordinates of the object by searching the provided name through the RFC catalog.
        """
        rfc_files = tuple((r.name for r in resources.files("vlbiplanobs.data").iterdir()
                           if r.is_file() and 'rfc' in r.name))
        assert len(rfc_files) > 0, "No RFC files found under the 'data' folder."

        with resources.as_file(resources.files("vlbiplanobs.data").joinpath(sorted(rfc_files)[-1])) \
                as rfcfile:
            process = subprocess.run(["grep", src_name, rfcfile], capture_output=True, text=True)
        # process = subprocess.run(["grep", src_name, "./data/rfc_2021c_cat.txt"], capture_output=True, text=True)

        if process.returncode == 1:
            raise ValueError(f"The source {src_name} was not found in the RFC catalog.")
        elif process.returncode != 0:
            raise RuntimeError("A problem happened while parsing the RFC catalog file.")
        elif process.stdout.count('\n') != 1:
            raise ValueError(f"None or multiple coincidences for the source {src_name} in the RFC catalog.")

        temp = process.stdout.split()  # Fields:  IVS  name  J2000name  h m s d m s
        return coord.SkyCoord(f"{temp[3]}h{temp[4]}m{temp[5]}s {temp[6]}d{temp[7]}m{temp[8]}s")

    def sun_separation(self, times: Time) -> Sequence[u.Quantity]:
        """Returns the separation of the source to the Sun at the given epoch(s).

        Inputs
        - times : astropy.time.Time
            An array of times defining the duration of the observation. That is,
            the first time will define the start of the observation and the last one
            will be the end of the observation. Note that the higher time resolution
            in `times` the more precise values you will obtain for antenna source
            visibility or determination of the rms noise levels. However, that will
            also imply a longer computing time.

        Return
        - astropy.units.Quantity
            The separation between the source and the Sun at the given times.
        """
        return self.coord.transform_to(coord.GCRS(obstime=times)).separation(coord.get_sun(times))

    def sun_constraint(self, min_separation: u.Quantity, times: Optional[Time] = None) -> Time:
        """Checks if the Sun can be a restriction to observe the given source.
        This is defined as if the Sun is at a separation lower than 'min_separation' from the
        target in the sky at a given moment.

        Inputs
        - min_separation : astropy.units.Quantity
            Minimun separation allowed for the Sun to be near the source.
            One can make use of the default separations provided in `freqsetups.solar_separations`
            for the different observing bands.
        - times : astropy.time.Time  [OPTIONAL]
            An array of times to check the separation of the Sun. If not provided, it will compute
            the full (current) calendar year with a resolution of one day.

        Returns
        - astropy.time.Time
            A list of times when the Sun is closer than the minimum separation.
            It would be an empty list if this condition never meets.
        """
        if times is None:
            times = Time(f"{dt.datetime.now(tz=dt.timezone.utc).year}-01-01 00:00") + \
                    np.arange(0, 365, 1)*u.day

        sun_separation = self.sun_separation(times=times)
        if isinstance(sun_separation, list):
            return [times[sun_separation[i] < min_separation] for i in range(len(self.coord))]
        else:
            return times[sun_separation < min_separation]


class SourceCatalog:
    def __init__(self, personal_catalog: Optional[str] = None):
        self._blocks: dict[str, dict[str, ScanBlock]] = dict()
        # self._sources: dict[str, Source] = dict()
        # self._all_names: dict[str, str] = dict()
        # self._personal: set[str] = set()
        # self._catalog: set[str] = set()
        if personal_catalog is not None:
            self.read_personal_catalog(personal_catalog)

    @property
    def blocknames(self):
        return [bb for b in self._blocks.values() for bb in b.keys()]

    @property
    def blocks(self):
        return {bb_key: bb_value for b in self._blocks.values() for bb_key, bb_value in b.items()}

    @property
    def targets(self):
        return self._blocks['targets']

    @property
    def pulsars(self):
        return self._blocks['pulsars'] if 'pulsars' in self._blocks else None

    @property
    def ampcals(self):
        return self._blocks['ampcals'] if 'ampcals' in self._blocks else None

    @property
    def fringefinders(self):
        return self._blocks['fringefinders'] if 'fringefinders' in self._blocks else None

    @property
    def polcals(self):
        return self._blocks['polcals'] if 'polcals' in self._blocks else None

    @functools.cache
    def source_names(self, include_calibrators: bool = False) -> list[str]:
        """Returns the names of all sources in the database.
        """
        if include_calibrators:
            return list([s.name for b in self._blocks.values() for bs in b.values() for s in bs.sources()])

        return list([s.name for b in self._blocks['targets'].values() for s in b.sources()])

    @functools.cache
    def sources(self, include_calibrators: bool = False) -> dict[str, Source]:
        """Returns all sources.
        """
        if include_calibrators:
            return {s.name: s for b in self._blocks.values() for bs in b.values() for s in bs.sources()}

        return {s.name: s for b in self._blocks['targets'].values() for s in b.sources()}

    def __contains__(self, item: str):
        return item in self.blocknames

    def __getitem__(self, item: str):
        for key in self._blocks:
            if item in self._blocks[key]:
                return self._blocks[key][item]

        raise KeyError(f"The item {item} is not in the catalog.")

    def read_personal_catalog(self, path: str):
        """Reads the yaml file with the information of all sources that may be scheduled.
        It returns the catalog as a dictionary.

        Input
            path : str
                The path to the yaml file with the catalog of sources to be imported.
        """
        with open(path, 'r') as sources_yaml:
            catalog = yaml.safe_load(sources_yaml)
            for a_entry in catalog:
                if a_entry not in self._blocks:
                    self._blocks[a_entry] = dict()

                match a_entry:
                    case 'pulsars':
                        src_type = SourceType.PULSAR
                    case 'targets':
                        src_type = SourceType.TARGET
                    case 'ampcals':
                        src_type = SourceType.AMPLITUDECAL
                    # case 'checksources':
                    #     src_type = SourceType.CHECKSOURCE
                    # case 'phasecals':
                    #     src_type = SourceType.PHASECAL
                    case 'fringefinders':
                        src_type = SourceType.FRINGEFINDER
                    case 'polcals':
                        src_type = SourceType.POLCAL
                    case _:
                        src_type = SourceType.UNKNOWN

                for a_src in catalog[a_entry]:
                    src = catalog[a_entry][a_src]
                    scans = []
                    if 'phasecal' in src:
                        scans.append(Scan(source=Source(name=src['phasecal']['name'],
                                                        coordinates=src['phasecal']['coordinates'],
                                                        source_type=SourceType.PHASECAL),
                                          duration=float(src['phasecal']['duration'])*u.min
                                          if 'duration' in src['phasecal'] else None,
                                          every=int(src['phasecal']['every'])
                                          if 'every' in src['phasecal'] else -1))

                    if 'checkSource' in src:
                        scans.append(Scan(source=Source(name=src['checkSource']['name'],
                                                        coordinates=src['checkSource']['coordinates'],
                                                        source_type=SourceType.CHECKSOURCE),
                                          duration=float(src['checkSource']['duration'])*u.min
                                          if 'duration' in src['checkSource'] else None,
                                          every=int(src['checkSource']['every'])
                                          if 'every' in src['checkSource'] else -1))

                    scans.append(Scan(source=Source(name=src['name'] if 'name' in src else a_src,
                                                    coordinates=src['coordinates'],
                                                    source_type=src_type),
                                      duration=float(src['duration'])*u.min if 'duration' in src else None,
                                      every=int(src['every']) if 'every' in src else -1))

                    self._blocks[a_entry][a_src] = ScanBlock(scans)

    def read_rfc_catalog(self, path: Optional[Union[str, Path]] = None):
        """Reads the catalog yaml file with the information of all sources that may be scheduled.
        It returns the catalog as a dictionary.

        Input
            path : str
                The path to the yaml file with the catalog of sources to be imported.
        """
        # TODO: convert this to another module and use duckDB, should be much faster
        raise NotImplementedError
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
    source: Source
    duration: u.Quantity = 10*u.min
    every: int = -1


class ScanBlock:
    """Defines a list of scans, each of them defined as a pointing to a given source during a given time.
    """

    def __init__(self, scans: list[Scan]):
        """Creates a block of scans, ideally a block to be observed with the target scans, phase-reference
        calibrator source (if needed), and check sources.

        For example, a block can consist on a single target scan (if this is a non phase-referencing
        observation).
        Then such scan will be repeated until fill the maximum observing time.

        On the contrary, in a phase-referencing observation a scan block can consist on one or
        multiple target scans, the associated phase-referencing scans, and possible check sources
        to be observed every certain number of target scans.

        If you want to schedule a fringe finder regularly during the observation, that is also possible.
        However, for isolated scans on these sources it will be preferred to be defined as different
        scan blocks.
        """
        if len(scans) == 0:
            raise ValueError("The scan block cannot contain an empty list of scans.")

        if not all([isinstance(s, Scan) for s in scans]):
            raise ValueError("All elements in the list of scans must be of Scan type.")

        self._scans = scans
        self._frac_time: dict[str, float] = {}

    @property
    def scans(self) -> list[Scan]:
        return self._scans

    def has(self, source_type: SourceType) -> bool:
        """Returns if the given source type is included among the ones observed in the provided
        list of scans.
        """
        return any([s.source.type is source_type for s in self.scans])

    def sources(self, source_type: Optional[SourceType] = None) -> list[Source]:
        """Returns the sources with the given source types in this block.
        """
        if source_type is None:
            return [s.source for s in self.scans]

        return [s.source for s in self.scans if s.source.type is source_type]

    def sourcenames(self, source_type: Optional[SourceType] = None) -> list[str]:
        """Returns the source names with the given source types in this block.
        """
        if source_type is None:
            return [s.source.name for s in self.scans]

        return [s.source.name for s in self.scans if s.source.type is source_type]

    def scan_with_sourcename(self, source_name: str) -> Optional[Scan]:
        for scan in self._scans:
            if scan.source.name == source_name:
                return scan

        raise ValueError(f"The source {source_name} is not present in any scan.")

    def scans_with_sources(self, source_type: SourceType) -> list[Scan]:
        """Returns the scans with the given source types in this block.
        """
        return [s for s in self._scans if s.source.type == source_type]

    def fractional_time(self) -> dict[str, float]:
        """Returns the fractional time dedicated to each source observed within the scan block.

        Returns
            fractional_time : dict[str, float]
                The fraction of time estimated to be spent on each particular source assuming
                the durations of the other scans to be observed in this block.
                The keys of the dict are the source names, and the values are the fraction of time,
                from the total scan block time, spent on the source.
        """
        if len(self._frac_time.keys()) > 0:
            return self._frac_time

        total_duration = sum([s.duration for s in self.scans if s.every <= 0])
        mcm_every = np.lcm.reduce([s.every for s in self.scans if s.every > 0])
        total_duration = total_duration*mcm_every + \
            sum([s.duration*(s.every/mcm_every) for s in self.scans if s.every > 0])

        for ascan in self.scans:
            self._frac_time[ascan.source.name] = ascan.duration*mcm_every / \
                                                 (ascan.every if ascan.every >= 0 else 1) / total_duration

        return self._frac_time

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

         NOTE: block scans should be easy!  This program is not prepared for the situation when you have
         multiple targets with multiple phase reference sources mixed.
        """
        # safety Checks
        if reduce(operator.add, [s.duration.to(u.min) for s in self.scans]) > max_duration:
            raise ValueError("The max_duration of the block cannot be shorter than the time of "
                             "all single scans.")

        main_loop: list[Scan] = []

        # First it arranges the (phasecal)/target scans if exist
        if self.has(SourceType.PHASECAL):
            if not self.has(SourceType.TARGET):
                raise ValueError("If phase calibrator scans provided, then target scans must also "
                                 "be provided.")

            loop_duration = reduce(operator.add,
                                   [s.duration for s in self.scans_with_sources(SourceType.TARGET) +
                                    self.scans_with_sources(SourceType.PHASECAL)])
            last_duration = reduce(operator.add,
                                   [s.duration for s in self.scans_with_sources(SourceType.PHASECAL)])
            # the second sum above is because the phase-referencing loop needs to be closed at the end.

            if any([s.every > -1 for s in self.scans
                    if s.source.type not in (SourceType.TARGET, SourceType.PHASECAL)]):
                # As there can be multiple sources to be observed every certain scans,
                # better to do it incremental...
                # full_n_reps = (max_duration - last_duration)/(loop_duration*)
                booked_time, n_loop = 0*u.min, 1
                while True:
                    target_in_this_scan: list[Scan] = []
                    for a_scan in [s for s in self.scans if s.every > -1]:
                        if n_loop % a_scan.every == 0:
                            target_in_this_scan += [a_scan,]

                    if len(target_in_this_scan) == 0:
                        target_in_this_scan = self.scans_with_sources(SourceType.TARGET)

                    to_append = self.scans_with_sources(SourceType.PHASECAL) + target_in_this_scan
                    n_loop += 1
                    if reduce(operator.add, [a.duration for a in to_append]) + \
                       booked_time > max_duration - last_duration:
                        break

                    main_loop += to_append
                    booked_time += reduce(operator.add, [a.duration for a in to_append])
            else:
                main_loop += (self.scans_with_sources(SourceType.PHASECAL) +
                              self.scans_with_sources(SourceType.TARGET)) * \
                              int(max_duration.to(u.min).value // (loop_duration +
                                  last_duration).to(u.min).value)

            main_loop += self.scans_with_sources(SourceType.PHASECAL)
        else:
            target_duration = reduce(operator.add,
                                     [s.duration for s in self.scans_with_sources(SourceType.TARGET)])
            main_loop += self.scans_with_sources(SourceType.TARGET) * \
                int(max_duration.to(u.min).value // target_duration.to(u.min).value)

        return main_loop

    def __iter__(self):
        yield from self._scans

    def __contains__(self, a_source_name: str):
        return a_source_name in self.sourcenames()
