import asyncio
from typing import Optional, Union, Iterable, Sequence, Self, Tuple
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
from .stations import Station, Stations
from .sources import Source, Scan, ScanBlock, SourceType, SourceNotVisible
from . import freqsetups


_NETWORKS = Stations.get_networks_from_configfile()
_STATIONS = Stations()


class Polarization(Enum):
    """Number of polarizations to be recorded. Possible values are:
    - SINGLE (only one polarization is recorded per antenna/baseline).
    - DUAL (the two RR,LL polarizations are recorded per antenna/baseline).
    - FULL (the full four RR,LL,RL,LR polarizations are recorded per antenna/baseline).
    """
    SINGLE = 1
    DUAL = 2
    FULL = 4

    def __str__(self):
        return self.name.lower()


class Observation(object):
    """Defines an observation of a single target source with a given network of stations.
    The observation can be set up in incremental steps (i.e. not all inputs are necessary
    at initialization time). Depending on the functions that you want to execute not all
    attributes may be required.
    An Observation allows the user to determine the elevation of the target source for
    each antenna in the network, the estimated rms noise level that would reach in such
    observation, or the resolution of the resulting images, assuming a standard neutral
    robust weighting.
    """

    def __init__(self, band: str, stations: Stations, scans: dict[str, ScanBlock],
                 times: Optional[Time] = None, duration: Optional[u.Quantity] = None,
                 datarate: Optional[Union[str, u.Quantity]] = None,
                 subbands: Optional[int] = None, channels: int = 64,
                 polarizations: Union[int, str] = 4,
                 inttime: Union[float, u.Quantity] = 2.0,
                 ontarget: float = 0.7, bits: int = 2):
        """Initializes an observation.
        Note that you can initialize an empty observation at this stage and add the
        information for the different attributes later. However, you may raise exception
        if running methods that require some of the unset attributes.

        Inputs
        - band : str
            Observing band to conduct the observation. Note that while this is a
            free-format string, it should match the type used when defining the bands
            at which each station can observe, and the default is the str format `XXcm`
            where XX represents the wavelength to observe in cm units.
            You can always check the available bands in `{Stations object}.observing_bands`.
        - stations : vlbiplanobs.stations.Stations
            Network of stations that will participate in the given observation.
        - scans : dict[str, ScanBlock]
            Scans to be observed during the observation.
            The key of the dict is just the name to refer to the specified ScanBlock as value.
        - times : astropy.time.Time
            An array of times defining the duration of the observation. That is,
            the first time will define the start of the observation and the last one
            will be the end of the observation. Note that the higher time resolution
            in `times` the more precise values you will obtain for antenna source
            visibility or determination of the rms noise levels. However, that will
            also imply a longer computing time. Steps of ~5-15 min seem appropriate
            for typical VLBI observations.
        - duration : astropy.units.Quantity
            Total duration of the observation. If 'times' provided, it will automatically be calculated
            from them. Otherwise it will be set from this value (e.g. for observations where the start
            time is not known).
        - datarate : int or astropy.units.Quantity
            Data rate for each antenna. It assumes that all antennas will run at the same
            data rate, which may not be true. If an int is introduce, it will be assumed to
            be given in Mbit/s. Otherwise, you can specify a Quantity with compatible units.
        - subbands : int
            Number of subbands in which the total bandwidth of the observation will be divided
            during correlation.
        - channels : int
            Number of channels for each subband to be created during correlation.
        - polarizations : int (1, 2, 4)  or str (single, dual, full)
            Number of polarizations to record in the observation. Three values are only possible
            0: single polarization.
            2: dual polarization (only direct terms: RR and LL, or XX and YY, are kept).
            4: full polarization (crossed terms are kept: RR, LL, RL, LR; or XX, YY, XY, YX).
        - inttime : float/int or astropy.units.Quantity
            Integration time (time resolution) that the final, correlated, data will show.
            If no units provided, seconds are assumed.
        - ontarget : float
            Fraction of the total observing time (end time minus start time) that will be
            spent on the target source. Note that in a typical VLBI observation only a fraction
            of the total observing time will end up in on-target source, commonly between 0.4-0.8
            (~40-80% of the time). This has an effect on the determination of the final rms noise level.
        - bits : int
            Number of bits at which the data have been recorded (sampled). A typical VLBI observation is
            almost always recorded with 2-bit data.
        """
        if isinstance(scans, dict) and all(isinstance(a_value, ScanBlock) for a_value in scans.values()):
            self.scans = scans
        else:
            raise ValueError("Scans needs to be either a ScanBlock or a list of ScanBlocks.")

        self._source_list: Optional[list[Source]] = None
        # Because otherwise gstimes is not initialized
        if times is None:
            self._times = None
            self._gstimes = None
        else:
            self.times = times

        self.band = band
        self.subbands = subbands
        self.channels = channels
        self.polarizations = polarizations if polarizations is not None \
            else Polarization.FULL  # type: ignore
        self.inttime = inttime
        self._duration = duration
        self.stations = stations
        self.bitsampling = bits
        self.datarate = datarate if isinstance(datarate, u.Quantity) or datarate is None \
            else datarate*u.Mbit/u.s
        self.ontarget_fraction = ontarget
        self._uv_baseline: Optional[dict[str, dict[str, u.Quantity]]] = None
        self._uv_array: Optional[dict[str, np.ndarray]] = None
        self._rms: Optional[Union[u.Quantity, dict[str, u.Quantity]]] = None
        self._synth_beam: Optional[dict[str, dict[str, u.Quantity]]] = None
        self._is_visible: Optional[dict[str, dict[str, coord.SkyCoord]]] = None
        self._is_always_visible: Optional[dict] = None
        self._altaz: Optional[dict] = None
        self._elevations: Optional[dict] = None

    @property
    def scans(self) -> Optional[dict[str, ScanBlock]]:
        """Returns the ScanBlock that will be observed during the observation.
        """
        return self._scans

    @scans.setter
    def scans(self, scans: Optional[dict[str, ScanBlock]]):
        self._scans = scans

    def sources(self, source_type: Optional[SourceType] = None) -> list[Source]:
        """Returns the sources included during the observation.
        It can return None if the target has not been set yet, showing a warning.

        Inputs
            source_type  : SourceType  (default = None)
                If set, it will only return the sources of type source_type.
                Otherwise, it returns all existing sources.

        Returns
            list[Source]
                List containing all existing sources matching the set criteria.
        """
        if self._scans is None:
            return []

        if self._source_list is not None:
            return self._source_list

        self._source_list = []
        for a_scanblock in self._scans.values():
            for a_src in a_scanblock.sources():
                if (source_type is None) or (a_src.type == source_type):
                    self._source_list.append(a_src)

        return self._source_list

    def _scanblock_name_from_source_name(self, source_name: str) -> Optional[str]:
        if self.scans is not None:
            for scanblockname, scanblock in self.scans.items():
                if source_name in scanblock:
                    return scanblockname
        else:
            return None

        raise ValueError(f"The source {source_name} is not present in any scan.")

    @property
    def sourcenames(self) -> list[str]:
        """Returns the source names included during the observation
        It can return None if the target has not been set yet, showing a warning.
        """
        return [a_source.name for a_source in self.sources()]

    @property
    def times(self) -> Optional[Time]:
        """Returns the times when the observation runs as an astropy.time.Time object.
        It can return None if the times have not been set yet, showing a warning.
        """
        return self._times

    @times.setter
    def times(self, new_times: Optional[Time]):
        if (new_times is not None) and (not isinstance(new_times, Time)):
            raise ValueError("'times' must be an astropy.time.Time instance or be None")
        elif isinstance(new_times, Time) and new_times.size < 2:
            raise ValueError("'times' must have at least two time values: "
                             "start and end of the observation.")

        self._times = new_times
        if self._times is None:
            self._gstimes = None
        else:
            self._gstimes = self._times.sidereal_time('mean', 'greenwich')

        self._uv_baseline = None
        self._uv_array = None
        self._rms = None
        self._synth_beam = None

    @property
    def gstimes(self) -> Optional[coord.angles.Longitude]:
        """Returns the GST times when the observation runs as an astropy.coordinates.angles.Longitude
        object (meaning in hourangle units).
        It can return None if the times have not been set yet, showing a warning.
        """
        return self._gstimes

    @property
    def duration(self) -> Optional[u.Quantity]:
        """Returns the total duration of the observation.
        """
        return self._duration if self.times is None else (self.times[-1] - self.times[0]).to(u.h)

    @property
    def band(self) -> str:
        """Returns the observing band at which the observation will be conducted.
        It can return None if the band has not been set yet, showing a warning.
        """
        return self._band

    @band.setter
    def band(self, new_band: str):
        if (new_band not in freqsetups.bands):
            raise ValueError("'new_band' needs to  match the following bands: "
                             f"{', '.join(freqsetups.bands.keys())} (wavelengths given in cm.)")

        self._band = new_band
        self._uv_baseline = None
        self._uv_array = None
        self._rms = None
        self._synth_beam = None

    @property
    def wavelength(self) -> u.Quantity:
        """Returns the central wavelength of the observation.
        """
        return float(self.band.replace('cm', ''))*u.cm

    @property
    def frequency(self) -> u.Quantity:
        """Returns the central frequency of the observations.
        """
        return 30*u.GHz/self.wavelength.to(u.cm).value

    @property
    def datarate(self) -> Optional[u.Quantity]:
        """Retuns the data rate (per station) used at which the observation is conducted.
        It can return None if the data rate has not been set yet, showing a warning.
        """
        return self._datarate

    @datarate.setter
    def datarate(self, new_datarate: Optional[u.Quantity]):
        """Sets the data rate used at each station during the observation.

        Inputs
        - new_datarate : astropy.units.Quantity  [e.g. Mb/s]
        """
        if new_datarate is None:
            self._datarate = None
        elif isinstance(new_datarate, int):
            if new_datarate <= 0:
                raise ValueError(f"datarate must be a positive number (currently {new_datarate})")

            self._datarate = new_datarate*u.Mbit/u.s
        elif isinstance(new_datarate, u.Quantity):
            if new_datarate <= 0:
                raise ValueError(f"datarate must be a positive number (currently {new_datarate})")

            self._datarate = new_datarate.to(u.Mbit/u.s)
        else:
            raise ValueError(f"Unknown type for datarate {new_datarate}"
                             "(int or astropy.units.Quantity (~bit/s) expected)")

        the_networks = self._guess_network()
        if 'EVN' in the_networks:
            for s in self.stations:
                if isinstance(s.max_datarate, dict) and 'EVN' in s.max_datarate:
                    s.datarate = min([d for d in (self._datarate, s.max_datarate['EVN'],
                                     _NETWORKS['EVN'].max_datarate(self.band)) if d is not None])
                else:
                    s.datarate = min([d for d in (self._datarate, s.max_datarate,
                                     _NETWORKS['EVN'].max_datarate(self.band)) if d is not None])
        else:
            for s in self.stations:
                s.datarate = min(self._datarate, s.max_datarate)

        if self._datarate is None:
            self._datarate = max([s.datarate for s in self.stations])

        if self.subbands is None:
            if 'EVN' in the_networks:  # 32-MHz subbands by default
                self.subbands = int(self._datarate.to(u.Mbit/u.s).value/32/8)
            else:  # 64-MHz subbannds by default
                self.subbands = int(self._datarate.to(u.Mbit/u.s).value/64/8)

        self._rms = None

    def _guess_network(self) -> list[str]:
        rating: dict = {}
        for network in _NETWORKS:
            if self.band in _NETWORKS[network].observing_bands:
                rating[network] = len([s for s in self.stations
                                       if s.codename in _NETWORKS[network].station_codenames]) / \
                                  len(self.stations)

        sorted_networks: list[str] = [n[0] for n in sorted(rating.items(), key=lambda x: rating[x[0]],
                                                           reverse=True)
                                      if rating[n[0]] > 0.0]
        return sorted_networks

    @property
    def subbands(self) -> Optional[int]:
        """Returns the number of subbands (also known as intermediate frequencies for AIPS users)
        in which the total bandwidth of the observation will be divided during correlation.
        It can return None if the number of subbands has not been set yet, showing a warning.
        """
        return self._subbands

    @subbands.setter
    def subbands(self, n_subbands: Optional[int]):
        if (n_subbands is not None) and (not isinstance(n_subbands, int)) and (n_subbands <= 0):
            raise ValueError("n_subbands needs to be a positive int, or None.")

        self._subbands = n_subbands

    @property
    def channels(self) -> int:
        """Returns the number of channels in which each subband will be divided during correlation.
        It can return None if the number of channels has not been set yet, showing a warning.
        """
        return self._channels

    @channels.setter
    def channels(self, n_channels: int):
        if n_channels is None:
            raise ValueError("channels needs to be a positive int, not None.")

        if (not isinstance(n_channels, int)) or (isinstance(n_channels, int) and n_channels <= 0):
            raise ValueError("channels needs to be a positive integer.")

        self._channels = n_channels

    @property
    def polarizations(self) -> Polarization:
        """Returns the number of polarizations that will be stored in the final data.
        It can return None if the number of polarizations has not been set yet, showing a warning.
        """
        return self._polarizations

    @polarizations.setter
    def polarizations(self, pols: Union[int, str, Polarization]):
        if isinstance(pols, Polarization):
            self._polarizations = pols
            return

        match pols:
            case 1 | 'single':
                self._polarizations = Polarization.SINGLE
            case 2 | 'dual':
                self._polarizations = Polarization.DUAL
            case 4 | 'full':
                self._polarizations = Polarization.FULL
            case _:
                raise ValueError("Polarizations needs to be either a Polarization value, "
                                 "the number 1, 2, or 4, or the strings "
                                 "'single', 'dual' or 'full' (equivalent to the numbers, respectively).")

    @property
    def inttime(self) -> Optional[u.Quantity]:
        """Returns the integration time used when correlating the observation as an astropy.units.Quantity.
        It can return None if the integration time has not been set yet, showing a warning.
        """
        return self._inttime

    @inttime.setter
    def inttime(self, new_inttime: Optional[Union[int, float, u.Quantity]]):
        """Sets the integration time of the observation.
        Inputs
        - new_inttime float/int or astropy.units.Quantity.
            If no units provided, seconds are assumed.
        """
        if new_inttime is None:
            self._inttime = None
        elif isinstance(new_inttime, float) or isinstance(new_inttime, int):
            if new_inttime <= 0:
                raise ValueError(f"'inttime' must be a positive number (currently {new_inttime})")

            self._inttime = new_inttime*u.s
        elif isinstance(new_inttime, u.Quantity):
            if new_inttime <= 0:
                raise ValueError(f"'inttime' must be a positive number (currently {new_inttime})")

            self._inttime = new_inttime.to(u.s)
        else:
            raise ValueError(f"Unknown type for 'inttime' {new_inttime} "
                             "(float/int/Quantity(~seconds) expected)")

    @property
    def ontarget_fraction(self) -> float:
        """Fraction of the total observing time spent on the target source.
        If scans include more sources than just the target source, then the fraction will be calculated
        from these scans, ignoring this fraction.
        """
        return self._ontarget

    @ontarget_fraction.setter
    def ontarget_fraction(self, ontarget: float):
        if not (0.0 < ontarget <= 1.0):
            raise ValueError("'ontarget_fraction' must be a float within (0.0, 1.0].")
        self._ontarget = ontarget
        self._rms = None

    @property
    def ontarget_time(self) -> Optional[dict[str, u.Quantity]]:
        """Total time spent on the target source during the observation.
        It can return None if the ontarget_fraction and duration have not been set yet, showing a warning.
        """
        if (self.duration is None) and (self.scans is None):
            return None
        elif (self.duration is not None) and (self.scans is None):
            return self.duration
        elif (self.duration is None) and (self.scans is not None):
            rprint("[red bold]Unexpected State[/red bold]: [red]You should have not reached "
                   "this place... Duration should be set already when scans are set![/red]")
            raise RuntimeError("Observation duration is not set but scans are set.")

        to_return = {}
        for scanblock in self.scans.values():  # type:ignore
            if len(scanblock.scans) == 1:
                to_return[scanblock.scans[0].source.name] = self.duration*self.ontarget_fraction  # type:ignore
            else:
                for src, dur in scanblock.fractional_time().items():
                    to_return[src] = dur * self.duration  # type: ignore

        return to_return

    @property
    def bandwidth(self) -> Optional[u.Quantity]:
        """Returns the total bandwidth of the observation.
        It returns None if the attributes 'polarizations', 'datarate', or 'bitsampling'
        have not been set yet.
        """
        if None in (self.polarizations, self.datarate, self.bitsampling):
            return None

        pols = self.polarizations.value % 3 + self.polarizations.value // 3  # Either 1 or 2
        return (self.datarate/(pols*self.bitsampling*2)).to(u.MHz)  # type: ignore

    @property
    def bitsampling(self) -> int:
        """Returns the bit sampling at which the data are recorded during the observation.
        """
        return self._bitsampling

    @bitsampling.setter
    def bitsampling(self, new_bitsampling: int):
        """Sets the bit sampling of the observation.
        Inputs
        - new_bitsampling : int
            In bits.
        """
        if isinstance(new_bitsampling, int):
            self._bitsampling = new_bitsampling*u.bit
        else:
            raise ValueError(f"Unknown type for {new_bitsampling} (int, in bits, expected)")

    @property
    def stations(self) -> Stations:
        """Returns the network of stations 'Stations' that will participate in this observation
        observing the target source.
        """
        return self._stations

    @stations.setter
    def stations(self, new_stations: Stations):
        assert isinstance(new_stations, Stations)
        self._stations = new_stations
        self._uv_baseline = None
        self._uv_array = None
        self._rms = None
        self._synth_beam = None
        self._elevations = None
        self._is_visible = None
        self._is_always_visible = None
        self._altaz = None

    def sources_in_block(self, block_name: str, source_type: Optional[SourceType] = None) -> list[Source]:
        """Returns the sources that are included in a specific scan block.
        """
        if self._scans is None or block_name not in self._scans:
            raise ValueError(f"The ScanBlock named {block_name} could not be found.")

        return self._scans[block_name].sources(source_type)

    def sourcenames_in_block(self, block_name: str, source_type: Optional[SourceType] = None) -> list[str]:
        """Returns the sources that are included in a specific scan block.
        """
        if self._scans is None or block_name not in self._scans:
            raise ValueError(f"The ScanBlock named {block_name} could not be found.")

        return self._scans[block_name].sourcenames(source_type) if self._scans is not None else None

    def elevations(self) -> dict[str, dict[str, dict[str, coord.angles.Latitude]]]:
        """Returns the elevation of the sources for each stations participating in the observation
        for all the given observing times.

        Returns
            elevations : dict
                Dictionary where they keys are: the ScanBlock name, the source name,
                and then the station code names, and the values are the elevations at each time stamp.
        """
        if (self.sources is None) or (self.times is None):
            raise ValueError("The target and/or observing times have not been initialized")
        # the current timing checks provided (on a test observation):
        # Elevations with process: 16.56423762498889 s
        # Elevations with threads: 0.2133850830141455 s
        # Elevations with nothing: 0.36229008401278406 s
        if self._elevations is None:
            self._elevations = self._elevations_threads()

        return self._elevations

    # INFO: this code in the following was meant to test how the computation runs faster
    # (process, threads, etc)
    def _elevation_func(self, a_scanblock: ScanBlock) -> dict[str, dict[str, Iterable]]:
        if (self.sources is None) or (self.times is None):
            raise ValueError("The target and/or observing times have not been initialized")

        # TODO: check if doing a multi-thread here will improve performance or not

        return {src.name: {a_station.codename: a_station.elevation(self.times, src)
                for a_station in self.stations}
                for src in self.sources()}

    def _elevations_process(self) -> dict:
        """Returns the elevation of the target source for each stations participating in the observation
        for all the given observing times.

        It computes it using different processes in parallel.

        Returns
            elevations : dict
                Dictionary where they keys are the station code names, and the values will be
                an astropy.coordinates.angles.Latitude object with the elevation at each observing time.
        """
        with ProcessPoolExecutor() as pool:
            results = pool.map(self._elevation_func, self.stations)

        return {s.codename: r for s, r in zip(self.stations, results)}

    def _elevations_threads(self) -> dict:
        """Returns the elevation of the target source for each stations participating in the observation
        for all the given observing times.

        It computes it using different concurrent threads.

        Returns
            elevations : dict
                Dictionary where they keys are the station code names, and the values will be
                an astropy.coordinates.angles.Latitude object with the elevation at each observing time.
        """
        assert self.scans is not None
        with ThreadPoolExecutor() as pool:
            results = pool.map(self._elevation_func, self.scans.values())

        return {s: r for s, r in zip(self.scans.keys(), results)}

    def _elevations_single(self) -> dict:
        """Returns the elevation of the target source for each stations participating in the observation
        for all the given observing times.

        It computes it without parallelization.

        Returns
            elevations : dict
                Dictionary where they keys are the station code names, and the values will be
                an astropy.coordinates.angles.Latitude object with the elevation at each observing time.
        """
        elevations = {}
        for a_station in self.stations:
            elevations[a_station.codename] = a_station.elevation(self.times, self.sources)
        return elevations

    async def _elevations_async_func(self, a_station: Station) -> list:
        return a_station.elevation(self.times, self.sources)

    async def _elevations_async_launch(self) -> dict:
        tasks: dict = {}
        for station in self.stations:
            tasks[station.codename] = asyncio.create_task(self._elevations_async_func(station))
            await tasks[station.codename]

        return tasks

    def _elevations_asyncio(self) -> dict:
        return asyncio.run(self._elevations_async_launch())
        # return {s.codename: v
        #         for s,v in zip(self.stations, asyncio.run(self._elevations_async_launch()))}

    def altaz(self) -> dict[str, dict[str, coord.SkyCoord]]:
        """Returns the altitude/azimuth of the target source for each stations participating
        in the observation for all the given observing times.

        Returns
            altaz : dict
                Dictionary where they keys are the station code names, and the values will be
                another dictionary with the source name as keys, and
                an astropy.coordinates.SkyCoord object with the altitude and azimuth
                of the source at each observing time.
        """
        if None in (self.sources, self.times):
            raise ValueError("The target and/or observing times have not been initialized")

        if self._altaz is None:
            self._altaz = self._altaz_threads()

        return self._altaz

    def _altaz_threads(self) -> dict:
        assert self.sources is not None and self.scans is not None

        def compute_altaz_for_source(station_source: tuple) -> tuple:
            station, source = station_source[0], station_source[1]
            return station, source, station.altaz(self.times, source)

        temp: dict[str, dict[str, coord.SkyCoord]] = {station.codename: {} for station in self.stations}
        for ablockname, ablock in self.scans.items():
            temp[ablockname] = {source.name: {} for source in ablock.sources()}
            with ThreadPoolExecutor() as executor:
                results = list(executor.map(compute_altaz_for_source, [(station, source)
                                            for station in self.stations for source in ablock.sources()]))
            for (station, source, altaz) in results:
                temp[ablockname][source.name][station.codename] = altaz

        return temp

    def _compute_is_visible_for_source(self, station_source_time) -> tuple:
        station, source = station_source_time[0], station_source_time[1]
        if len(station_source_time) == 3:
            times = station_source_time[2]
        else:
            times = self.times
        return station, source, station.is_observable(times, source)

    def is_observable(self) -> dict[str, dict[str, list[bool]]]:
        """Returns whenever the given ScanBlock can be observed by each station for each time
        of the observation.

        Returns
            is_visible : dict
                Dictionary where they keys are the station code names, and the values will be
                a dictionary with the sources names as keys, and a tuple containing
                a numpy array with the indexes in the `Observation.times`
                array with the times where the target source can be observed by the station.

                In this sense, you can e.g. call obs.times[obs.is_visible[a_station_codename]]
                to get such times.
        """
        if None in (self.scans, self.times, self.sources):
            raise ValueError("The target, sources and/or observing times have not been initialized")

        if self._is_visible is None:
            self._is_visible = {}
            for ablockname, ablock in self.scans.items():  # type: ignore
                with ThreadPoolExecutor() as executor:
                    results = list(executor.map(self._compute_is_visible_for_source, [(station, source)
                                                for station in self.stations
                                                for source in ablock.sources()]))  # type: ignore

                self._is_visible[ablockname] = {}
                for (station, source, visible) in results:
                    if station.codename in self._is_visible[ablockname]:
                        self._is_visible[ablockname][station.codename] *= visible
                    else:
                        self._is_visible[ablockname][station.codename] = visible

        return self._is_visible

    def is_observable_at(self, time: Time) -> dict[str, dict[str, list[bool]]]:
        """Returns whenever the given ScanBlock can be observed by each station at the given time.
        """
        _is_visible_at: dict[str, dict[str, list[bool]]] = {}
        for ablockname, ablock in self.scans.items():  # type: ignore
            with ThreadPoolExecutor() as executor:
                results = list(executor.map(self._compute_is_visible_for_source, [(station, source, time)
                                            for station in self.stations
                                            for source in ablock.sources()]))  # type: ignore

            _is_visible_at[ablockname] = {}
            for (station, source, visible) in results:
                _is_visible_at[ablockname][station.codename] = visible

        return _is_visible_at

    def can_be_observed(self) -> dict[str, dict[str, bool]]:
        """Returns whenever the sources can be observed by each station at least during part
        of the observation.

        Returns
            is_observable : dict
                Dictionary where they keys are the station code names, and the values will be
                a dictionary with the sources names as keys, and a boolean containing
                if the source can be observed by the station.
        """
        return {s: {k: any(v) for k, v in i.items()} for s, i in self.is_observable().items()}

    def is_always_observable(self) -> dict[str, dict[str, bool]]:
        """Returns whenever the target source can be observed by each station along the whole observation.

        Returns
            is_visible : dict
                Dictionary where they keys are the station code names, and the values will be
                a tuple containing a numpy array with the boolean indicating if the given BlockScan
                (same order as in Observation.sources) is always visible for the given station.
        """
        if None in (self.scans, self.times, self.sources):
            raise ValueError("The target, sources and/or observing times have not been initialized")

        if self._is_always_visible is None:
            def compute_is_always_visible_for_source(station_source: tuple) -> tuple:
                station, source = station_source[0], station_source[1]
                return station, source, station.is_always_observable(self.times, source)

            self._is_always_visible = {ablockname: {} for ablockname in self.scans}  # type: ignore
            for ablockname, ablock in self.scans.items():  # type: ignore
                with ThreadPoolExecutor() as executor:
                    results = list(executor.map(compute_is_always_visible_for_source, [(station, source)
                                                for station in self.stations
                                                for source in ablock.sources()]))  # type: ignore

                for (station, source, visible) in results:
                    self._is_always_visible[ablockname][station.codename] = visible

        return self._is_always_visible

    def when_is_observable(self, min_stations: int = 3,
                           mandatory_stations: Optional[str | list[str]] = None,
                           stations_all_time: bool = False, within_time_range: Optional[Time] = None,
                           return_gst: bool = False) -> dict[str, list[Time | u.Quantity]]:
        """Returns the time range when the different BlockScans/sources are visible verifying the set
        requirements.
        If the time of the observation is set, it will check if the source is observable only within that
        time range, otherwise it will check at any possible time.

        Inputs
        - min_stations : int  [default = 3]
            Minimum number of stations that are allowed to consider the source as observable.
        - mandatory_stations : list[str]  [default = None]
            Defines the stations that must be able to observe the source to be consider for observations.
            It must be a list of strings with the codes or names of the stations.
            The wildcard 'all' is possible, indicating that only times where all antennas can observe
            the source are allowed.
        - stations_all_time : bool  [default = False]
            If all stations must be observing the source at all times or this is not required.
        - within_time_range : Time  [default = None]
            Time range within which the source is observable.
            Must contain only two times (start and end time).
        - return_gst : bool  [default = False]
            Defines if the returned times for the start and end of the observation should be in GST time
            instead of UTC.

        Returns
            - (t0, t1, ...)  : tuple with two astropy.units.Quantity[Longitude] objects.
                The start and end (GST) time of the same period of time when the source is visible by
                enough stations.
                If 'return_gst' is False, the tuple will contain two astropy.time.Time objects instead,
                with the UTC times for the start and end times.
                It will be always an even number of entries. If more than two, means that there are
                multiple time ranges where it can be observed following t0, t1, t2, t3, ...
                so the source can be observed between t0 and t1, t2 and t3, etc.

        Exceptions
        - It may raise the exception SourceNotVisible if the target source is not visible by
          enough stations.
        """
        if self.sources is None:
            raise ValueError("The sources have not been initialized")

        if within_time_range is None:
            if self.times is None:
                # Using the autonm equinox as that one GST ~= UTC
                within_time_range = Time('2025-09-21', scale='utc') + np.arange(0.0, 1.005, 0.01)*u.day
            else:
                within_time_range = self.times
        else:
            assert len(within_time_range) == 2, \
                "The time range must contain only two times (start and end time)"
            assert within_time_range[1] <= within_time_range[0], \
                "The end time must be larger than the start time"
            within_time_range = within_time_range[0] + \
                np.arange((within_time_range[1] - within_time_range[0]).to(u.day).value, 0.01)*u.day

        if mandatory_stations == 'all':
            mandatory_stations = self.stations.station_codenames
            min_stations = len(mandatory_stations)

        result: dict[str, list[Time | u.Quantity]] = {}
        for blockname, station_visibility in self.is_observable_at(time=within_time_range).items():  # type: ignore
            visible_times = [False for i in range(len(within_time_range))]
            for i, time in enumerate(within_time_range):
                if ((sum(vis[i] for vis in station_visibility.values()) >= min_stations)
                    and ((mandatory_stations is None)
                         or all(station_visibility[station][i] for station in mandatory_stations))):
                    visible_times[i] = True

            diff = [i + 1*(not value) for i, value in enumerate(visible_times)
                    if value != (visible_times + [False])[i+1] or (i == 0 and value)]
            result[blockname] = list(zip(within_time_range[diff[::2]], within_time_range[diff[1::2]]))
            # It may happen that the first and last times are contiguous
            if len(result[blockname]) > 1 and \
               np.abs(result[blockname][0][0].mjd - result[blockname][-1][1].mjd) % 1 < 0.01:
                result[blockname][-1] = (result[blockname][-1][0], result[blockname][0][1] + 1*u.day)
                result[blockname] = result[blockname][1:]

        if return_gst:
            return {s: list([tuple(t.sidereal_time('mean', 'greenwich') for t in tu)
                    for tu in result[s]]) for s in result}

        return result

    def scheduler(self) -> list[Scan]:
        """For the given observing time and scan blocks, it will schedule in an optimal way the scans
        for the different sources along the observation.
        """
        raise NotImplementedError

    def longest_baseline(self) -> dict[str, Tuple[str, u.Quantity]]:
        """Returns the longest baseline in the observation for each source.

        Returns
        - 'source name': ('{ant1}-{ant2}', length) : dict[str, tuple[str, astropy.units.Quantity]]
            - '{ant1}-{ant2}' : str
                Composed by the codenames of the two antennas (ant1, ant2) conforming the longest baseline.
            - length : astropy.units.Quantity
                The projected length of the baseline as seen from the target source position.
        """
        uv_data = self.get_uv_data()
        longest_bl: dict[str, Tuple[str, u.Quantity]] = {}
        for src, src_uv in uv_data.items():
            max_temp = -1.0
            for a_bl, uv in src_uv.items():
                if len(uv) > 0:
                    bl_length = np.sqrt(np.max((uv**2).sum(axis=1)))
                    if bl_length > max_temp:
                        longest_bl[src] = (a_bl, (bl_length*self.wavelength).to(u.km))
                        max_temp = bl_length

        return longest_bl

    def shortest_baseline(self) -> dict[str, Tuple[str, u.Quantity]]:
        """Returns the shortest baseline in the observation.

        Returns
        - 'source name': ('{ant1}-{ant2}', length) : tuple
            - '{ant1}-{ant2}' : str
                Composed by the codenames of the two antennas (ant1, ant2) conforming the
                shortest baseline.
            - length : astropy.units.Quantity
                The projected length of the baseline as seen from the target source position.
        """
        uv_data = self.get_uv_data()
        shortest_bl: dict[str, Tuple[str, u.Quantity]] = {}
        for src, src_uv in uv_data.items():
            min_temp = -1
            for a_bl, uv in src_uv.items():
                if len(uv) > 0:
                    bl_length = np.sqrt(np.min((uv**2).sum(axis=1)))
                    if bl_length < min_temp or min_temp < 0:
                        shortest_bl[src] = (a_bl, (bl_length*self.wavelength).to(u.km))
                        min_temp = bl_length

        return shortest_bl

    def bandwidth_smearing(self) -> u.Quantity:
        """Returns the bandwidth smearing expected for the given observation.

        The peak response to a point target source decreases at positions farther away from the
        pointing (correlated) sky position due to the frequency averaging performed in the data.

        This function returns the angular separation at which the bandwidth smearing produces
        a reduction of a 10% in the response of the telescope. The field of view should then
        be limited to this range to avoid significant loses.
        """
        return ((49500*u.arcsec*u.MHz*u.km)*self.channels /
                (np.max([lb[1] for lb in self.longest_baseline().values()]) *
                 self.bandwidth/self.subbands)).to(u.arcsec)

    def time_smearing(self) -> u.Quantity:
        """Returns the time smearing expected for the given observation.

        The peak response to a point target source decreases at positions farther away from the
        pointing (correlated) sky position due to the time averaging performed in the data.

        This function returns the angular separation at which the time smearing produces
        a reduction of a 10% in the response of the telescope. The field of view should then
        be limited to this range to avoid significant loses.
        """
        return ((18560*u.arcsec*u.km*u.s/u.cm) *
                (self.wavelength / (np.max([lb[1]
                 for lb in self.longest_baseline()])*self.inttime))).to(u.arcsec)

    def datasize(self) -> Optional[u.Quantity]:
        """Returns the expected size for the output FITS IDI files.

        A regular observation with the European VLBI Network is stored in FITS IDI files,
        typically several 2-GB files. This function provides an estimation of the total
        size for these stored files.
        Note that this function does not take into account down times for the different
        stations. The provided value will thus always be un upper-limit for the real, final,
        value.
        """
        if None in (self.duration, self.inttime, self.subbands, self.channels):
            return None

        temp = len(self.stations)**2 * (self.duration/self.inttime).decompose()  # type: ignore
        temp *= self.polarizations.value*self.subbands*self.channels  # type: ignore
        return temp*1.75*u.GB/(131072*3600)

    def thermal_noise(self) -> Optional[Union[u.Quantity, dict[str, u.Quantity]]]:
        """Returns the expected rms thermal noise for the given observation.

        Each antenna has a different sensitivity for the observing band (established from
        the SEFD value). This function computes the available baselines at each timestamp
        and estimates the final thermal noise reached when integrating the whole observation.

        Note that this function takes into account when each station can observe the source,
        but does not take into account sensitivity drops doe to external factors like e.g.
        low elevations of the source. The provided thermal noise is also assumed when a natural
        weighting is applied to the data when imaging. The thermal noise can thus be a bit
        higher if an uniform roubst is used.

        If the source has not been set, it will assume that all antennas observe all time.
        """
        if self._rms is not None:
            return self._rms

        if self.datarate is None:
            raise ValueError("The data rate must be defined to estimate thermal noise.")

        sefds = np.array([stat.sefd(self.band).to(u.Jy).value for stat in self.stations])
        bandwidths = np.array([s.datarate.to(u.bit/u.s).value if s.datarate is not None
                              else self.datarate.to(u.bit/u.s).value for s in self.stations]) * \
            2 * 2 * min(self.polarizations.value, 2)  # In Hz
        bandwidth_min = np.minimum.outer(bandwidths, bandwidths)
        if len(self.sources()) == 0 or self.times is None:
            if self.times is not None:
                dt = (self.times[-1]-self.times[0]).to(u.min).value
            elif self.duration is not None:
                dt = self.duration.to(u.min).value
            else:
                raise ValueError("Either observinng times or duration must be defined "
                                 "to estimate thermal noise.")

            # THIS IS THE NEW IMPLEMENTATION TO USE.
            temp = np.sum((bandwidth_min/np.outer(sefds, sefds)) *
                          np.triu(np.ones_like(bandwidth_min), k=1))

            # THIS IS THE FASTEST IMPLEMENTATION BUT NEEDS TO INCLUDE THE DIFFERENT BANDWIDTHS
            # return (np.sum(1/(sefds[:, None]*sefds[None, :])) - np.sum(1/(sefds*sefds))) / 2
            # temp = (np.sum(((bandwidth_min/np.outer(sefds, sefds)) *
            #                 np.triu(np.ones_like(bandwidth_min), k=1)) *
            #                (1/(sefds[:, None]*sefds[None, :]) - 1/(sefds*sefds)))) / 2
            # OLD IMPLEMENTATION
            # temp = 0.0
            # for j in range(len(sefds)):
            #     for k in range(j+1, len(sefds)):
            #         temp += dt/(sefds[j]*sefds[k])

            # Quicker:  sefds = np.array(..)
            # return (np.sum(1/sefds[:, None] * 1/sefds[None, :]) - np.sum(1/sefds ** 2)) / 2
            # This runs much faster: 46 us versus 24 ms.

            # temp = 1.0/np.sqrt(temp*self.ontarget_fraction)
            self._rms = ((np.sqrt(2)/0.7)/np.sqrt(temp*dt))*u.Jy/u.beam
            return self._rms
        else:
            self._rms = {}
            delta_t = (self.times[1] - self.times[0]).to(u.min).value
            for sourcename in self.sourcenames:
                scanname = self._scanblock_name_from_source_name(sourcename)
                assert scanname is not None, f"No scan found related to the source {sourcename}"
                visible = np.array([self.is_observable()[scanname][stat.codename] for stat in self.stations])
                integrated_time = np.sum(visible[:, None] & visible[None, :], axis=2) * delta_t
                temp = np.sum((bandwidth_min*integrated_time /
                              np.outer(sefds, sefds)) * np.triu(np.ones_like(bandwidth_min), k=1))
                self._rms[sourcename] = ((np.sqrt(2)/0.7)/np.sqrt(temp))*u.Jy/u.beam

            return self._rms

    def _compute_uv_per_source(self, source: Optional[Source] = None) -> dict[str, u.Quantity]:
        nstat = len(self.stations)
        # Determines the xyz of all baselines. Time independent
        bl_xyz = np.empty(((nstat*(nstat-1))//2, 3))
        bl_names = []
        s = [ant.location for ant in self.stations]
        for i in range(nstat):
            for j in range(i+1, nstat):
                # An unique number defining a baseline
                k = int(i*(nstat-1) - sum(range(i)) + j-i)
                bl_xyz[k-1, :] = np.array([ii.value for ii in s[i].to_geocentric()]) - \
                    np.array([ii.value for ii in s[j].to_geocentric()])
                bl_names.append("{}-{}".format(self.stations[i].codename,
                                               self.stations[j].codename))

        if source is None:
            hourangle: coord.Angle = self.gstimes
            print("WARNING: 'target' is not set, thus we assume a source at +/- 45º declination"
                  " to estimate the (u, v) values.'")
            m = np.array([[np.sin(hourangle), np.cos(hourangle), np.zeros(len(hourangle))],
                          [-np.sin(45*u.deg)*np.cos(hourangle),
                           np.sin(45*u.deg)*np.sin(hourangle),
                           np.cos(45*u.deg)*np.ones(len(hourangle))]])
        else:
            hourangle = (self.gstimes - source.ra.to(u.hourangle)).value % 24*u.hourangle
            m = np.array([[np.sin(hourangle), np.cos(hourangle), np.zeros(len(hourangle))],
                          [-np.sin(source.dec)*np.cos(hourangle),
                           np.sin(source.dec)*np.sin(hourangle),
                           np.cos(source.dec)*np.ones(len(hourangle))]])

        bl_uv = np.array([m[:, :, i] @ bl_xyz.T for i in range(m.shape[-1])])*u.m
        bl_uv_up: dict[str, u.Quantity] = {}
        if source is None:
            for i, bl_name in enumerate(bl_names):
                ant1, ant2 = bl_name.split('-')
                bl_uv_up[bl_name] = (bl_uv[:, :, i] / self.wavelength).decompose()

            return bl_uv_up
        else:
            ants_up = self.is_observable()
            blockname = [ablock for ablock in self.scans
                         if source.name in self.scans[ablock].sourcenames()][0]
            for i, bl_name in enumerate(bl_names):
                ant1, ant2 = bl_name.split('-')
                bl_up = (np.array([a for a in ants_up[blockname][ant1]
                                   if a in ants_up[blockname][ant2]]), )
                if len(bl_up[0]) > 0:
                    bl_uv_up[bl_name] = (bl_uv[:, :, i][bl_up] / self.wavelength).decompose()

            if len(bl_uv_up.keys()) == 0:
                raise SourceNotVisible

        return bl_uv_up

    def get_uv_data(self) -> dict[str, dict[str, u.Quantity]]:
        """Returns the (u, v) values for each baseline and each timestamp for which the source
        is visible.

        Returns
        - {source name: '{ant1}-{ant2}': uv_data} : dict
            - '{ant1}-{ant2}' : str
                The keys of the dictionary are the baselines forming the array, where ant1 and ant2
                are the code names of each station.
            - uv_data : astropy.units.Quantity
                A 2-d array with all (u, v) values for each timestamp in `times` when the source is
                visible for the given baseline. (u,v) are given in lambda units.
                Note that complex conjugate values are not provided.

        Exceptions
        - It may raise the exception SourceNotVisible if no baselines can observe the source
          at all during the observation.
        """
        if self._uv_baseline is not None:
            return self._uv_baseline

        with ThreadPoolExecutor() as executor:
            self._uv_baseline = dict(zip([src.name for src in self.sources()],
                                         executor.map(self._compute_uv_per_source, self.sources())))

        return self._uv_baseline

    def get_uv_values(self) -> dict[str, np.ndarray]:
        """Returns the (u, v) values for each baseline and each timestamp for which the source
        is visible.

        The difference with `get_uv_baseline` is that `get_uv_array` only returns the (u,v)
        values, dropping the information of baselines and times to which these values belong to.

        Returns
            uv_values : dict[str, np.ndarray]
            The dict keys are the source names, and the values a (N, 2)-dimensional
            numpy.ndarray containing all N (u,v) points resulting for
            each timestamp and baseline. The (u,v) values are given in lambda units.
            Note that complex conjugate values are not provided.

        Exceptions
        - It may raise the exception SourceNotVisible if no baselines can observe the source
          at all during the observation.
        """
        if self._uv_array is not None:
            return self._uv_array

        self._uv_array = {}
        bl_uv_up = self.get_uv_data()
        for src, bl_uv in bl_uv_up.items():
            self._uv_array[src] = np.empty((np.sum([bl_uv[bl_name].shape[0] for bl_name in bl_uv]), 2))
            last_i = 0
            for bl_name in bl_uv:
                self._uv_array[src][last_i:last_i+bl_uv[bl_name].shape[0], :] = bl_uv[bl_name]
                last_i += bl_uv[bl_name].shape[0]

        return self._uv_array

    def synthesized_beam(self) -> dict[str, dict[str, u.Quantity]]:
        """Estimates the resulting synthesized beam of the observations based on
        the expected (u,v) coverage per source.

        This is just an estimation made by a ellipse fitting to the (u, v) coverage,
        from which we obtain the resolution on the two axes following
        https://science.nrao.edu/facilities/vlba/docs/manuals/oss/ang-res
            theta_HPBW (mas)  \\sim 2063 x lambda(cm)/b_max^km

        Note that different robust weighting maps during imaging would provide slightly
        different synthesized beams. The provided value here does not assumed any weighting
        in the data. A natural weighting would thus likely provide a slightly larger beam,
        while an uniform weighting would provide a slightly smaller beam.

        Returns a dict with the following keys: 'bmaj', 'bmin', 'pa' (major axis, minor axis,
        and position angle).
        The three values are astropy.units.Quantity objects with units of angles.
        """
        if self._synth_beam is not None:
            return self._synth_beam

        self._synth_beam = {}

        resolution = lambda bl : ((2.063e8*u.mas)/bl).to(u.mas)
        uvvis = self.get_uv_values()
        for src, uv in uvvis.items():
            # Transform the uv points into r,theta (polar) points
            uvvis_polar = np.empty_like(uv)
            uvvis_polar[:,0] = np.sqrt((uv**2).sum(axis=1)) # radius
            uvvis_polar[:,1] = np.arctan2(uv[:,1], uv[:,0]) # theta
            # Defines the BMAJ and PA
            bl_bmaj = np.max(uvvis_polar[:,0])
            bl_bmaj_theta = uvvis_polar[:,1][np.where(uvvis_polar[:,0] == bl_bmaj)][0]
            # Gets the BMIN and an orthogonal projection
            bl_bmin_theta = ( bl_bmaj_theta + np.pi/2 ) % (2*np.pi)
            bl_bmin = np.max(np.abs(uv.dot(np.array([np.cos(bl_bmin_theta),
                                                     np.sin(bl_bmin_theta)]))))

            self._synth_beam[src] = {'bmaj': resolution(bl_bmin), 'bmin': resolution(bl_bmaj),
                                      'pa': (bl_bmaj_theta*u.rad).to(u.deg)}

        return self._synth_beam


    def get_dirtymap(self, pixsize: int = 1024, robust: str = "natural", oversampling: int = 4):
        """Returns the dirty beam produced for the given observation.

        Input:
            - pixsize : int
                Size in pixels of the returned (squared) image. By default 1024.
            - robust : str
                The weighting Briggs robust used to compute the dirty map.
                It must be either 'natural' or 'uniform'.
            - oversampling : int
                Oversampling factor when plotting the dirty map.
                Recommended values are between 1 and 10.
                CURRENTLY THIS PARAMETER IS IGNORED. WILL BE IMPLEMENTED IN A LATER VERSION

        Returns:
            - dirty_image : (pixsize x pixsize) np.array
                The dirty image in intensity.
            - laxis : np.array
                An array representing the values of each pixel in the image in mas for each axis.

        """
        assert robust in ('natural', 'uniform')
        assert isinstance(pixsize, int)
        # assert isinstance(oversampling, int) and 20 > oversampling >= 1
        # Oversampling the dirty beam
        uvimg = np.zeros((pixsize, pixsize))

        # Gridding uv.  Inspired from VNSIM code (https://github.com/ZhenZHAO/VNSIM)
        # uvdata = self.get_uv_array()
        # uvscaling = pixsize/(2.0*1.025*np.max(np.max(np.abs(uvdata), axis=0)))
        # uvimg[pixsize//2 + np.trunc(uvdata[:,0]*uvscaling).astype(int), \
        #       pixsize//2 + np.trunc(uvdata[:,1]*uvscaling).astype(int)] += 1
        #
        # # Recovering the requested size for the image
        # # if oversampling > 1:
        # #     uvimg = scipy.ndimage.zoom(uvimg, oversampling, order=1)
        # # uvimg = uvimg[int(pixsize*0.5):int(pixsize*1.5), int(pixsize*0.5):int(pixsize*1.5)]
        # if robust == 'uniform':
        #     if np.max(uvimg) > 1:
        #         uvimg[np.where(uvimg > 0)] = 1
        #     else:
        #         print('\n\nWARNING: Uniform dirty map same as natural one.\n\n')
        #
        # dirty_beam = np.real(np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(uvimg/np.max(uvimg)))))
        # if oversampling > 1:
        #     dirty_beam = scipy.ndimage.zoom(dirty_beam, oversampling, order=1)
        # imgsize = (uvscaling*u.rad/(2*oversampling)).to(u.mas) # angular equivalent size of the resulting image
        # return dirty_beam[int(pixsize*(oversampling-1)/2):int(pixsize*(oversampling+1)/2), \
        #                   int(pixsize*(oversampling-1)/2):int(pixsize*(oversampling+1)/2)].T/np.max(dirty_beam), \
        #        np.linspace(-imgsize, imgsize, pixsize)
        #         # return dirty_beam.T/np.max(dirty_beam), \
        #

        oversampling = 1
        # Gridding uv.  Inspired from VNSIM code (https://github.com/ZhenZHAO/VNSIM)
        uvdata = self.get_uv_array()
        uvscaling = (oversampling*pixsize)/(2*1.05*np.max(np.max(np.abs(uvdata), axis=0)))
        if robust == 'natural':
            uvimg[oversampling*pixsize//2 + np.trunc(uvdata[:,0]*uvscaling).astype(int), \
                  oversampling*pixsize//2 + np.trunc(uvdata[:,1]*uvscaling).astype(int)] += 1
        else:
            uvimg[oversampling*pixsize//2 + np.trunc(uvdata[:,0]*uvscaling).astype(int), \
                  oversampling*pixsize//2 + np.trunc(uvdata[:,1]*uvscaling).astype(int)] = 1

        # Recovering the requested size for the image
        # uvimg = uvimg[int(pixsize*0.5):int(pixsize*1.5), int(pixsize*0.5):int(pixsize*1.5)]
        dirty_beam = np.real(np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(uvimg/np.max(uvimg)))))
        imgsize = (uvscaling*u.rad/2).to(u.mas) # angular equivalent size of the resulting image
        return dirty_beam[int(pixsize*(oversampling-1)/2):int(pixsize*(oversampling+1)/2), \
                          int(pixsize*(oversampling-1)/2):int(pixsize*(oversampling+1)/2)].T/np.max(dirty_beam), \
               np.linspace(-imgsize, imgsize, pixsize)

    def print_obs_times(self, date_format='%d %B %Y'):
        """Returns the time range (starttime-to-endtime) of the observation in a smart way.
        If the observation lasts for less than one day it omits the end date:
                20 January 1971
                10:00-20:00 UTC
                GST range: 05:00-15:00

        If the observation ends the day after, then it returns:
                20 January 1971
                10:00-20:00 UTC (+1d)
                GST range: 05:00-15:00

        If the observation is longer, then it returns
                20 January 1971 10:00 to  24 January 1971 20:00 UTC
                GST range: 05:00-15:00

        Input:
            - date_format : str [OPTIONAL]
                Format for the date part (only the date part, not the times) of
                the string to represent the time range.
        Output:
            - strtime : str
                A string showing the time-range of the observation.
        """
        if self.gstimes is None:
            return "Observing time unspecifed."
        # TODO: this can be written better with self.gstimes[0].to_string(sep=':', precision=2, pad=True)
        gsttext = "{:02n}:{:02.2n}-{:02n}:{:02.2n}".format((self.gstimes[0].hour*60) // 60,
                                                           (self.gstimes[0].hour*60) % 60,
                                                           (self.gstimes[-1].hour*60) // 60,
                                                           (self.gstimes[0].hour*60) % 60)
        if self.times[0].datetime.date() == self.times[-1].datetime.date():
            return "{} {}-{} UTC\nGST range: {}".format(self.times[0].datetime.strftime(date_format),
                                                        self.times[0].datetime.strftime('%H:%M'),
                                                        self.times[-1].datetime.strftime('%H:%M'), gsttext)
        elif (self.times[-1] - self.times[0]) < 24*u.h:
            return "{} {}-{} UTC (+1d)\nGST range: {}".format(
                                        self.times[0].datetime.strftime(date_format),
                                        self.times[0].datetime.strftime('%H:%M'),
                                        self.times[-1].datetime.strftime('%H:%M'), gsttext)
        else:
            return "{} {} to {} {} UTC\nGST range: {}".format(
                                        self.times[0].datetime.strftime(date_format),
                                        self.times[0].datetime.strftime('%H:%M'),
                                        self.times[-1].datetime.strftime(date_format),
                                        self.times[-1].datetime.strftime('%H:%M'), gsttext)


    def get_fig_ant_elev(self):
        data_fig = []
        data_dict = self.elevations()
        # Some reference lines at low elevations
        for ant in data_dict:
            data_fig.append({'x': self.times.datetime, 'y': data_dict[ant].value,
                            'mode': 'lines', 'hovertemplate': "Elev: %{y:.2n}º<br>%{x}",
                            'name': self.stations[ant].name})

        data_fig.append({'x': self.times.datetime, 'y': np.zeros_like(self.times)+10,
                         'mode': 'lines', 'hoverinfo': 'skip', 'name': 'Elev. limit 10º',
                         'line': {'dash': 'dash', 'opacity': 0.5, 'color': 'gray'}})
        data_fig.append({'x': np.unwrap(self.gstimes.value*2*np.pi/24)*24/(2*np.pi), 'y': np.zeros_like(self.times)+20,
                         'xaxis': 'x2', 'mode': 'lines', 'hoverinfo': 'skip',
                         'name': 'Elev. limit 20º', 'line': {'dash': 'dot', 'opacity': 0.5,
                         'color': 'gray'}})
        return {'data': data_fig,
                'layout': {'title': 'Source elevation during the observation',
                           'hovermode': 'closest',
                           'xaxis': {'title': 'Time (UTC)', 'showgrid': False,
                                     'ticks': 'inside', 'showline': True, 'mirror': False,
                                     'hovermode': 'closest', 'color': 'black'},
                           'xaxis2': {'title': {'text': 'Time (GST)', 'standoff': 0},
                                      'showgrid': False, 'overlaying': 'x', #'dtick': 1.0,
                                      'tickvals': np.arange(np.ceil(self.gstimes.value[0]),
                                            np.floor(np.unwrap(self.gstimes.value*2*np.pi/24)[-1]*24/(2*np.pi))+1),
                                      'ticktext': np.arange(np.ceil(self.gstimes.value[0]),
                                            np.floor(np.unwrap(self.gstimes.value*2*np.pi/24)[-1]*24/(2*np.pi))+1)%24,
                                      'ticks': 'inside', 'showline': True, 'mirror': False,
                                      'hovermode': 'closest', 'color': 'black', 'side': 'top'},
                           'yaxis': {'title': 'Elevation (degrees)', 'range': [0., 92.],
                                     'ticks': 'inside', 'showline': True, 'mirror': "all",
                                     'showgrid': False, 'hovermode': 'closest'},
                           'zeroline': True, 'zerolinecolor': 'k'}}




    def get_fig_ant_up(self):
        data_fig = []
        data_dict = self.is_visible()
        gstimes = np.unwrap(self.gstimes.value*2*np.pi/24)*24/(2*np.pi)
        gstimes = np.array([dt.datetime(self.times.datetime[0].year, self.times.datetime[0].month,
                            self.times.datetime[0].day) + dt.timedelta(seconds=gst*3600) for gst in gstimes])
        for i,ant in enumerate(data_dict):
            # xs = [obs.times.datetime[0].date() + datetime.timedelta(seconds=i*3600) for i in np.unwrap(obs.gstimes.value*2*np.pi/24)[data_dict[ant]]*24/(2*np.pi)]
            xs = gstimes[data_dict[ant]]
            data_fig.append({'x': xs,
                             'y': np.zeros_like(data_dict[ant][0])-i, 'type': 'scatter',
                             'hovertemplate': "GST %{x}",
                             'mode': 'markers', 'marker_symbol': "41",
                             'hoverinfo': "skip",
                             'name': self.stations[ant].name})

        data_fig.append({'x': self.times.datetime, 'y': np.zeros_like(self.times)-0.5,
                         'xaxis': 'x2',
                         'mode': 'lines', 'hoverinfo': 'skip', 'showlegend': False,
                         'line': {'dash': 'dot', 'opacity': 0.0, 'color': 'white'}})
        return {'data': data_fig,
                'layout': {'title': {'text': 'Source visible during the observation',
                                     'y': 1, 'yanchor': 'top'},
                           'hovermode': 'closest',
                           'xaxis': {'title': 'Time (GST)', 'showgrid': False,
                                     'range': [gstimes[0], gstimes[-1]],
                                     # 'tickvals': np.arange(np.ceil(obs.gstimes.value[0]),
                                     #        np.floor(np.unwrap(obs.gstimes.value*2*np.pi/24)[-1]*24/(2*np.pi))+1),
                                     # 'ticktext': np.arange(np.ceil(obs.gstimes.value[0]),
                                     #        np.floor(np.unwrap(obs.gstimes.value*2*np.pi/24)[-1]*24/(2*np.pi))+1)%24,
                                     'tickformat': '%H:%M',
                                     'ticks': 'inside', 'showline': True, 'mirror': False,
                                     'hovermode': 'closest', 'color': 'black'},
                           'xaxis2': {'title': {'text': 'Time (UTC)', 'standoff': 0},
                                      'showgrid': False, 'overlaying': 'x', #'dtick': 1.0,
                                      'ticks': 'inside', 'showline': True, 'mirror': False,
                                      'hovermode': 'closest', 'color': 'black', 'side': 'top'},
                           'yaxis': {'ticks': '', 'showline': True, 'mirror': True,
                                     'showticklabels': False, 'zeroline': False,
                                     'showgrid': False, 'hovermode': 'closest',
                                     'startline': False}}}



    def get_fig_uvplane(self):
        data_fig = []
        bl_uv = self.get_uv_baseline()
        for bl_name in bl_uv:
            # accounting for complex conjugate
            uv = np.empty((2*len(bl_uv[bl_name]), 2))
            uv[:len(bl_uv[bl_name]), :] = bl_uv[bl_name]
            uv[len(bl_uv[bl_name]):, :] = -bl_uv[bl_name]
            data_fig.append({'x': uv[:,0],
                             'y': uv[:,1],
                             # 'type': 'scatter', 'mode': 'lines',
                             'type': 'scatter', 'mode': 'markers',
                             'marker': {'symbol': '.', 'size': 2},
                             'name': bl_name, 'hovertext': bl_name, 'hoverinfo': 'name', 'hovertemplate': ''})
        return {'data': data_fig,
                'layout': {'title': '', 'showlegend': False,
                           'hovermode': 'closest',
                           'width': 700, 'height': 700,
                           'xaxis': {'title': 'u (lambda)', 'showgrid': False, 'zeroline': False,
                                     'ticks': 'inside', 'showline': True, 'mirror': "all",
                                     'color': 'black'},
                           'yaxis': {'title': 'v (lambda)', 'showgrid': False, 'scaleanchor': 'x',
                                     'ticks': 'inside', 'showline': True, 'mirror': "all",
                                     'color': 'black', 'zeroline': False}}}



    def get_fig_dirty_map(self):
        # Right now I only leave the natural weighting map (the uniform does not always correspond to the true one)
        dirty_map_nat, laxis = self.get_dirtymap(pixsize=1024, robust='natural', oversampling=4)
        raise NotImplementedError
        # TODO: uncomment these two lines and move them outside observation. Flexibility
        # fig1 = px.imshow(img=dirty_map_nat, x=laxis, y=laxis[::-1], labels={'x': 'RA (mas)', 'y': 'Dec (mas)'}, \
        #         aspect='equal')
        # fig = make_subplots(rows=1, cols=1, subplot_titles=('Natural weighting',), shared_xaxes=True, shared_yaxes=True)
        fig.add_trace(fig1.data[0], row=1, col=1)
        mapsize = 30*self.synthesized_beam()['bmaj'].to(u.mas).value
        fig.update_layout(coloraxis={'showscale': False, 'colorscale': 'Inferno'}, showlegend=False,
                          xaxis={'autorange': False, 'range': [mapsize, -mapsize]},
                          yaxis={'autorange': False, 'range': [-mapsize, mapsize]}, autosize=False)
        fig.update_xaxes(title_text="RA (mas)", constrain="domain")
        fig.update_yaxes(title_text="Dec (mas)", scaleanchor="x", scaleratio=1)
        # dirty_map_nat, laxis = obs.get_dirtymap(pixsize=1024, robust='natural', oversampling=4)
        # dirty_map_uni, laxis = obs.get_dirtymap(pixsize=1024, robust='uniform', oversampling=4)
        # fig1 = px.imshow(img=dirty_map_nat, x=laxis, y=laxis[::-1], labels={'x': 'RA (mas)', 'y': 'Dec (mas)'}, \
        #         aspect='equal')
        # fig2 = px.imshow(img=dirty_map_uni, x=laxis, y=laxis[::-1], labels={'x': 'RA (mas)', 'y': 'Dec (mas)'}, \
        #         aspect='equal')
        # fig = make_subplots(rows=1, cols=2, subplot_titles=('Natural weighting', 'Uniform weighting'),
        #                     shared_xaxes=True, shared_yaxes=True)
        # fig.add_trace(fig1.data[0], row=1, col=1)
        # fig.add_trace(fig2.data[0], row=1, col=2)
        # mapsize = 30*obs.synthesized_beam()['bmaj'].to(u.mas).value
        # fig.update_layout(coloraxis={'showscale': False, 'colorscale': 'Inferno'}, showlegend=False,
        #                   xaxis={'autorange': False, 'range': [mapsize, -mapsize]},
        #                   # This xaxis2 represents the xaxis for fig2.
        #                   xaxis2={'autorange': False, 'range': [mapsize, -mapsize]},
        #                   yaxis={'autorange': False, 'range': [-mapsize, mapsize]}, autosize=False)
        # fig.update_xaxes(title_text="RA (mas)", constrain="domain")
        # fig.update_yaxes(title_text="Dec (mas)", row=1, col=1, scaleanchor="x", scaleratio=1)

        return fig

    # def export_to_pdf(self, outputname):
    #     """Exports the basic information of the observation into a PDF file.
    #     """
    # https://towardsdatascience.com/creating-pdf-files-with-python-ad3ccadfae0f
    #     from fpdf import FPDF
    #     pdf = FPDF(orientation='P', unit='cm', format='A4')
    #     # pdf_w = 21.0
    #     # pdf_h = 29.7
    #     pdf.add_page()
    #
    #     plotly.io.write_image(pltx,file='pltx.png',format='png',width=700, height=450)
    #     pltx=(os.getcwd()+'/'+"pltx.png")
    #
    #     pdf.image(sctplt,  link='', type='', w=1586/80, h=1920/80)
    #     self.set_font('Arial', 'B', 16)
    #     self.set_text_color(220, 50, 50)
    #     self.cell(w=210.0, h=40.0, align='C', txt="LORD OF THE PDFS", border=0)
    #     self.set_font('Arial', '', 12)
    #     self.multi_cell(0,10,txt)
    #
    #     pdf.output('test.pdf','F')
    #


