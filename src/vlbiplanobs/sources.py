from typing import Optional, Union, Self
from importlib import resources
import subprocess
# from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import numpy as np
import operator
from functools import reduce
import datetime as dt
from enum import Enum, auto
from dataclasses import dataclass
from astropy import units as u
from astropy.time import Time
from astropy import coordinates as coord
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


class SourceType(Enum):
    """Types of sources in a regular VLBI observation
    """
    TARGET = auto()
    PHASECAL = auto()
    FRINGEFINDER = auto()
    AMPLITUDECAL = auto()
    CHECKSOURCE = auto()
    POLCAL = auto()
    UNKNOWN = auto()


class Source(FixedTarget):
    """Defines a target source located at some coordinates and with a given name.
    """
    def __init__(self, name: Optional[str] = None,
                 coordinates: Optional[str] = None,
                 source_type: SourceType = SourceType.UNKNOWN, **kwargs):
        """Initializes a Source object.

        Inputs
        - name : str  [OPTIONAL]
            Name associated to the source. By default is None.
        - coordinates : str [OPTIONAL]
            Coordinates of the target source in a str format recognized by
            astropy.coordinates.SkyCoord (e.g. XXhXXmXXs XXdXXmXXs).
            J2000 coordinates are assumed.
            If not provided, name must be a source name recognized by the RFC catalog or astroquery.
        - source_type : SourceType  [default = UNKNOWN]
            Defines the type of the source.
        - kwargs
            keyword arguments to be passed to astropy.coordinates.SkyCoord() if needed.
            For example, the 'unit=' parameter.

        If both provided, the given coordinates will be used for the given source.

        It may raise:
        - NameResolveError: if the name is not recognized (and no coordinates are provided)
        - ValueError: if the coordinates have an unrecognized format.
        - AttributeError: if neither name or coordinates are provided, or name is empty.
        """
        if ((name is None) and (coordinates is None)) or (name == ''):
            raise AttributeError("At least one 'coordiantes' or 'name' (not empty) must be provided")

        if not isinstance(source_type, SourceType):
            raise ValueError("source_type must be a SourceType value.")

        if (name is not None) and (coordinates is None):
            coordinates = self.get_coordinates_from_name(name)

        super().__init__(coord.SkyCoord(coordinates, **kwargs), name)

        self._type = source_type

    @property
    def type(self) -> SourceType:
        return self._type

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
    def get_source_from_name(cls, src_name: str) -> Self:
        """Returns a Source object by finding the coordinates from its name.
        It first searches wihtin the RFC catalog, and if not found, through the
        'astropy.coordinates.get_icrs_coordinates' function through the RFC catalogs.
        """
        return cls(src_name, coordinates=Source.get_coordinates_from_name(src_name))


    @staticmethod
    def get_rfc_coordinates(src_name: str) -> coord.SkyCoord:
        """Returns the coordinates of the object by searching the provided name through the RFC catalog.
        """
        rfc_files = tuple((r.name for r in resources.files("data").iterdir() if r.is_file() and 'rfc' in r.name))
        assert len(rfc_files) > 0, "No RFC files found under the 'data' folder."

        with resources.as_file(resources.files("data").joinpath(sorted(rfc_files)[-1])) as rfcfile:
            process = subprocess.run(["grep", src_name, rfcfile], capture_output=True, text=True)

        if process.returncode == 1:
            raise ValueError(f"The source {src_name} was not found in the RFC catalog.")
        elif process.returncode != 0:
            raise RuntimeError("A problem happened while parsing the RFC catalog file.")
        elif process.stdout.count('\n') != 1:
            raise ValueError(f"None or multiple coincidences for the source {src_name} in the RFC catalog.")

        temp = process.stdout.split()  # Fields:  IVS  name  J2000name  h m s d m s
        return coord.SkyCoord(f"{temp[3]}h{temp[4]}m{temp[5]}s {temp[6]}d{temp[7]}m{temp[8]}s")


    @classmethod
    def get_rfc_source(cls, src_name: str) -> Self:
        """Returns a Source object by searching the provided name through the RFC catalog.
        """
        return cls(src_name, Source.get_rfc_coordinates(src_name))


    def sun_separation(self, times: Time) -> Union[list, np.array]:
        """Returns the separation of the source to the Sun at the given epoch(s).

        Inputs
        - times : astropy.time.Time
            An array of times defining the duration of the observation. That is,
            the first time will define the start of the observation and the last one
            will be the end of the observation. Note that the higher time resolution
            in `times` the more precise values you will obtain for antenna source
            visibility or determination of the rms noise levels. However, that will
            also imply a longer computing time.
        """
        if len(self.coord.shape) > 0:
            return [c.transform_to(coord.GCRS(obstime=times)).separation(coord.get_sun(times)) \
                    for c in self.coord]
        else:
            return self.coord.transform_to(coord.GCRS(obstime=times)).separation(coord.get_sun(times))


    def sun_constraint(self, min_separation: u.Quantity, times: Optional[Time] = None) -> Union[list, np.array]:
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


@dataclass
class Scan:
    source: Source
    duration: u.Quantity
    every: int = -1


class ScanBlock:
    """Defines a list of scans, each of them defined as a pointing to a given source during a given time.
    """
    def __init__(self, scans: list[Scan]):
                 #target: Union[tuple, Scan], phasecal: Optional[Union[tuple, Scan]] = None,
                 # checksource: Optional[Union[tuple, Scan]] = None,
                 # fringefinder: Optional[Union[tuple, Scan]] = None,
                 # duration: Optional[u.Quantity] = None,
                 # max_duration: Optional[u.Quantity] = None, min_duration: Optional[u.Quantity] = None,
                 # wavelength: Optional[u.Quantity] = None):
        """Creates a block of scans, ideally a block to be observed with the target scans, phase-reference
        calibrator source (if needed), and check sources.

        For example, a block can consist on a single target scan (if this is a non phase-referencing observation).
        Then such scan will be repeated until fill the maximum observing time.

        On the contrary, in a phase-referencing observation a scan block can consist on one or multiple target
        scans, the associated phase-referencing scans, and possible check sources to be observed every certain
        number of target scans.

        If you want to schedule a fringe finder regularly during the observation, that is also possible.
        However, for isolated scans on these sources it will be preferred to be defined as different scan blocks.
        """
        if len(scans) == 0:
            raise ValueError("The scan block cannot contain an empty list of scans.")

        if not all([isinstance(s, Scan) for s in scans]):
            raise ValueError("All elements in the list of scans must be of Scan type.")

        self._scans = scans


    @property
    def scans(self) -> list[Scan]:
        return self._scans


    def has(self, source_type: SourceType) -> bool:
        """Returns if the given source type is included among the ones observed in the provided list of scans.
        """
        return any([s.source.type is source_type for s in self.scans])

    def sources(self, source_type: Optional[SourceType] = None) -> list[Source]:
        """Returns the sources with the given source types in this block.
        """
        if source_type is None:
            return [s.source for s in self.scans]

        return [s.source for s in self.scans if s.source.type is source_type]


    def fill(self, max_duration: u.Quantity) -> list[Scan]:
        """Given the list of scans, returns the final arrangement of scans that fills the available time.
        This will follow the following conditions:
        - Repeats the target scan within the given time. If a `phasecal` is provided, then it will always
          bracket each target scan with this `phasecal`. If two or more phasecal are provided (e.g. P1, P2),
          then it will assume a multi phase-referencing technique for the target T:  P1 P2 T P1 P2 T P1 P2...
        - If some sources have the 'every = N > 0' constraint, then these will be scheduled every N cycles.
          For example, if no phasecal are provided and N = 3 for the C source, then: T T C T T C ...
          And if one phasecal is provided: P T P T P C P T P...

         NOTE: block scans should be easy!  This program is not prepared for the situation when you have
         multiple targets with multiple phase reference sources mixed.
        """
        # safety Checks
        if reduce(operator.add, [s.duration.to(u.min) for s in self.scans]) > max_duration:
            raise ValueError("The max_duration of the block cannot be shorter than the time of all single scans.")

        main_loop: list[Scan] = []


        # First it arranges the (phasecal)/target scans if exist
        if self.has(SourceType.PHASECAL):
            if not self.has(SourceType.TARGET):
                raise ValueError("If phase calibrator scans provided, then target scans must also be provided.")

            loop_duration = reduce(operator.add, [s.duration for s in self.sources(SourceType.TARGET) + \
                                                                      self.sources(SourceType.PHASECAL)])
            last_duration = reduce(operator.add, [s.duration for s in self.sources(SourceType.PHASECAL)])
            # the second sum above is because the phase-referencing loop needs to be closed at the end.

            if any([s.every > -1 for s in self.scans \
                                 if s.source.type not in (SourceType.TARGET, SourceType.PHASECAL)]):
                # As there can be multiple sources to be observed every certain scans, better to do it incremental...
                # full_n_reps = (max_duration - last_duration)/(loop_duration*)
                booked_time, n_loop = 0*u.min, 1
                while True:
                    target_in_this_scan: list[Source] = []
                    for a_scan in [s for s in self.scans if s.every > -1]:
                        if n_loop % a_scan.every == 0:
                            target_in_this_scan += [a_scan.source]

                    if len(target_in_this_scan) == 0:
                        target_in_this_scan = self.sources(SourceType.TARGET)

                    to_append = self.sources(SourceType.PHASECAL) + target_in_this_scan
                    n_loop += 1
                    if reduce(operator.add, [a.duration for a in to_append]) + booked_time > \
                                                                 max_duration - last_duration:
                        break

                    main_loop += to_append
                    booked_time += reduce(operator.add, [a.duration for a in to_append])
            else:
                main_loop += (self.sources(SourceType.PHASECAL) + self.sources(SourceType.TARGET)) * \
                             int(max_duration.to(u.min).value // (loop_duration + last_duration).to(u.min).value)

            main_loop += self.sources(SourceType.PHASECAL)
        else:
            target_duration = reduce(operator.add, [s.duration for s in self.sources(SourceType.TARGET)])
            main_loop += self.sources(SourceType.TARGET) * \
                         int(max_duration.to(u.min).value // target_duration.to(u.min).value)

        return main_loop






    # def __getitem__(self, item: Union[str, int]):
    #     """Returns the source associated to the given index position or source name.
    #     """
    #     if isinstance(item, str):
    #         return Source(name=item, coordinates=self.coord[self.name.index(item)])
    #     elif isinstance(item, int):
    #         return Source(name=self.name[item], coordinates=self.coord[item])
    #     else:
    #         raise ValueError("'item' needs to be either a string with the name of the source or " \
    #                 "the index in the current source list.")

    # def __len__(self):
    #     if len(self.coord.shape) > 0:
    #         return self.coord.shape[0]
    #     else:
    #         return 1


    # def __iter__(self):
    #     # TODO: I think when source == 1, then this goes into an infinite loop
    #     if len(self.coord.shape) > 0:
    #         yield Source(self.name[self._index], self.coord[self._index])
    #     else:
    #         yield Source(self.name, self.coord)


    # def __next__(self):
    #     self._index += 1
    #     if self._index > len(self):
    #         raise StopIteration
    #     # if self._index > len(self.coord.shape):
    #     # The standard action is to increase the current counter (e.g. self._index, = 0 created in init), and return the element at that index. If the index > len, then raise StopIteration.


