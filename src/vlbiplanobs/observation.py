import threading
from functools import wraps
from collections import defaultdict
from typing import Optional, Union, Tuple, Literal, get_type_hints
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from enum import Enum
from astropy import units as u
from astropy import coordinates as coord
from astropy.time import Time
from .stations import Stations, Station
from .sources import Source, Scan, ScanBlock, SourceType, SourceNotVisible
from . import freqsetups


_NETWORKS = Stations.get_networks_from_configfile()
_STATIONS = Stations()


def enforce_types(func):
    """Decorator that will raise TypeError if the passed attribute to a given function
    has the wrong type.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        hints = get_type_hints(func)
        # Check positional dash_bootstrap_components
        for arg, (name, expected_type) in zip(args, hints.items()):
            if name == 'return':
                continue

            if not isinstance(arg, expected_type):
                raise TypeError(f"Argument '{name}' must be of type {expected_type.__name__}, "
                                f"but got {type(arg).__name__}")

        # Check kwyword arguments
        for name, arg in kwargs.items():
            if name not in hints:
                continue

        expected_type = hints[name]
        if not isinstance(arg, expected_type):
            raise TypeError(f"Argument '{name}' must be of type {expected_type.__name__}, "
                            f"but got {type(arg).__name__}")

        return func(*args, **kwargs)

    return func



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

    @property
    def _REF_TIMES(self):
        return Time('2025-09-21', scale='utc') + np.arange(0.0, 1.005, 0.01)*u.day

    @property
    def _REF_YEAR(self):
        return Time('2025-01-01', scale='utc') + np.arange(0.0, 365.2, 1)*u.day

    @enforce_types
    def __init__(self, band: str, stations: Stations, scans: dict[str, ScanBlock],
                 times: Optional[Time] = None, duration: u.Quantity = 24*u.h,
                 datarate: Optional[u.Quantity] = None,
                 subbands: int = 8, channels: int = 64,
                 polarizations: Union[int, str] = 4,
                 inttime: u.Quantity = 2.0*u.s,
                 ontarget: float = 0.7, bits: u.Quantity = 2*u.bit):
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
            time is not known). If not given it will default to 24 h.
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
        self._mutex: threading.Lock = threading.Lock()
        self._mutex_uv: threading.Lock = threading.Lock()
        if isinstance(scans, dict) and all(isinstance(a_value, ScanBlock) for a_value in scans.values()):
            self.scans = scans
        else:
            raise ValueError("Scans needs to be either a ScanBlock or a list of ScanBlocks.")

        self._source_list: Optional[list[Source]] = None
        self._fixed_time: bool = times is not None
        self.times = times
        self.band = band
        self.subbands = subbands
        self.channels = channels
        self.polarizations = polarizations if polarizations is not None \
            else Polarization.FULL  # type: ignore
        self.inttime = inttime
        self._duration = duration if duration is not None else 24*u.h
        self.stations = stations
        self.bitsampling = bits
        self.datarate = datarate
        self.ontarget_fraction = ontarget
        self._uv_baseline: Optional[dict[str, dict[str, u.Quantity]]] = None
        self._uv_array: Optional[dict[str, np.ndarray]] = None
        self._rms: Optional[Union[u.Quantity, dict[str, u.Quantity]]] = None
        self._synth_beam: Optional[dict[str, dict[str, u.Quantity]]] = None
        self._is_visible: Optional[dict[str, dict[str, coord.SkyCoord]]] = None
        self._is_always_visible: Optional[dict] = None
        self._altaz: Optional[dict] = None
        self._elevations: Optional[dict] = None
        self._baseline_sensitivity: Optional[dict[str, u.Quantity]] = None

    @property
    def scans(self) -> dict[str, ScanBlock]:
        """Returns the ScanBlock that will be observed during the observation.
        """
        return self._scans

    @scans.setter
    @enforce_types
    def scans(self, scans: Optional[dict[str, ScanBlock]]):
        with self._mutex:
            self._scans = scans if scans is not None else {}
            self._elevations = None
            self._altaz = None
            self._is_visible = None
            self._is_always_observable = None
            self._rms = None
            self._uv_baseline = None
            self._uv_array = None
            self._synth_beam = None

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

    def _scanblock_name_from_source_name(self, source_name: str) -> str:
        if self.scans is not None:
            for scanblockname, scanblock in self.scans.items():
                if source_name in scanblock:
                    return scanblockname

        raise ValueError(f"The source {source_name} is not present in any scan.")

    @property
    def sourcenames(self) -> list[str]:
        """Returns the source names included during the observation
        It can return None if the target has not been set yet, showing a warning.
        """
        return [a_source.name for a_source in self.sources()]

    @property
    def times(self) -> Time:
        """Returns the times when the observation runs as an astropy.time.Time object.
        It can return None if the times have not been set yet, showing a warning.
        """
        return self._times

    @times.setter
    @enforce_types
    def times(self, new_times: Optional[Time]):
        if (new_times is not None) and isinstance(new_times, Time) and new_times.size < 2:
            raise ValueError("'times' must have at least two time values: "
                             "start and end of the observation.")
        with self._mutex:
            self._times = new_times if new_times is not None else self._REF_TIMES
            self._gstimes = self._times.sidereal_time('mean', 'greenwich')
            self._fixed_time = new_times is not None

            self._elevations = None
            self._altaz = None
            self._is_visible = None
            self._is_always_observable = None
            self._rms = None
            self._uv_baseline = None
            self._uv_array = None
            self._synth_beam = None

    @property
    def gstimes(self) -> coord.angles.Longitude:
        """Returns the GST times when the observation runs as an astropy.coordinates.angles.Longitude
        object (meaning in hourangle units).
        It can return None if the times have not been set yet, showing a warning.
        """
        return self._gstimes

    @property
    def fixed_time(self) -> bool:
        return self._fixed_time

    @property
    def duration(self) -> u.Quantity:
        """Returns the total duration of the observation.
        """
        return self._duration if not self.fixed_time else (self.times[-1] - self.times[0]).to(u.h)

    @duration.setter
    @enforce_types
    def duration(self, new_duration: u.Quantity):
        """Sets the total duration of the observation.
        """
        if not new_duration.unit.is_equivalent(u.h):
            raise TypeError("The new duration must be a quantity with time units.")

        self._duration = new_duration
        with self._mutex:
            self._rms = None
            self._baseline_sensitivity = None
            self._uv_baseline = None
            self._is_always_observable = None
            self._uv_array = None
            self._synth_beam = None

    @property
    def band(self) -> str:
        """Returns the observing band at which the observation will be conducted.
        It can return None if the band has not been set yet, showing a warning.
        """
        return self._band

    @band.setter
    @enforce_types
    def band(self, new_band: str):
        if (new_band not in freqsetups.bands):
            raise TypeError("'new_band' needs to  match the following bands: "
                            f"{', '.join(freqsetups.bands.keys())} (wavelengths given in cm.)")

        with self._mutex:
            self._band = new_band
            self._rms = None
            self._baseline_sensitivity = None
            self._uv_baseline = None
            self._uv_array = None
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
    def datarate(self) -> u.Quantity:
        """Retuns the data rate (per station) used at which the observation is conducted.
        It can return None if the data rate has not been set yet, showing a warning.
        """
        return self._datarate

    @datarate.setter
    @enforce_types
    def datarate(self, new_datarate: u.Quantity):
        """Sets the data rate used at each station during the observation.

        Inputs
        - new_datarate : astropy.units.Quantity  [e.g. Mb/s]
        """
        if not new_datarate.unit.is_equivalent(u.bit/u.s):
            raise TypeError("The new data rate must be a quantity with bit /s equivalent units.")

        if new_datarate <= 0:
            raise ValueError(f"datarate must be a positive number (currently {new_datarate})")

        the_networks = self._guess_network()
        with self._mutex:
            self._datarate = new_datarate.to(u.Mbit/u.s)
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
                    s.datarate = min(self._datarate, s.max_datarate) if s.max_datarate is not None \
                                 else self._datarate

            if self._datarate is None:
                self._datarate = max([s.datarate for s in self.stations])

            if self.subbands is None:
                if 'EVN' in the_networks:  # 32-MHz subbands by default
                    self.subbands = int(self._datarate.to(u.Mbit/u.s).value/32/8)
                else:  # 64-MHz subbannds by default
                    self.subbands = int(self._datarate.to(u.Mbit/u.s).value/64/8)

            self._rms = None
            self._baseline_sensitivity = None

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

    @staticmethod
    def guess_network(band: str, list_stations: list[Station]) -> list[str]:
        """Returns the VLBI network that is driving the observations, making a guess given
        the antennas that are participating and the observing band.
        """
        rating: dict = {}
        for network in _NETWORKS:
            if band in _NETWORKS[network].observing_bands:
                rating[network] = len([s for s in list_stations
                                       if s.codename in _NETWORKS[network].station_codenames]) / \
                                  len(list_stations)

        sorted_networks: list[str] = [n[0] for n in sorted(rating.items(), key=lambda x: rating[x[0]],
                                                           reverse=True)
                                      if rating[n[0]] > 0.0]
        return sorted_networks

    @staticmethod
    def _get_max_datarate(network: str, band: str) -> u.Quantity:
        return _NETWORKS[network].max_datarate(band)

    @property
    def subbands(self) -> int:
        """Returns the number of subbands (also known as intermediate frequencies for AIPS users)
        in which the total bandwidth of the observation will be divided during correlation.
        It can return None if the number of subbands has not been set yet, showing a warning.
        """
        return self._subbands

    @subbands.setter
    @enforce_types
    def subbands(self, n_subbands: int):
        if n_subbands <= 0:
            raise ValueError(f"n_subbands needs to be a positive integer, but is {n_subbands}.")

        self._subbands = n_subbands

    @property
    def channels(self) -> int:
        """Returns the number of channels in which each subband will be divided during correlation.
        It can return None if the number of channels has not been set yet, showing a warning.
        """
        return self._channels

    @channels.setter
    @enforce_types
    def channels(self, n_channels: int):
        if n_channels <= 0:
            raise ValueError(f"n_channels needs to be a positive integer, but is {n_channels}.")

        self._channels = n_channels

    @property
    def polarizations(self) -> Polarization:
        """Returns the number of polarizations that will be stored in the final data.
        It can return None if the number of polarizations has not been set yet, showing a warning.
        """
        return self._polarizations

    @polarizations.setter
    @enforce_types
    def polarizations(self, pols: Union[int, str, Polarization]):
        with self._mutex:
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
                                     "the number 1, 2, or 4, or the strings 'single', "
                                     "'dual' or 'full' (equivalent to the numbers, respectively).")

            self._rms = None
            self._baseline_sensitivity = None
            self._uv_baseline = None
            self._uv_array = None
            self._synth_beam = None

    @property
    def inttime(self) -> u.Quantity:
        """Returns the integration time used when correlating the observation as an astropy.units.Quantity.
        It can return None if the integration time has not been set yet, showing a warning.
        """
        return self._inttime

    @inttime.setter
    @enforce_types
    def inttime(self, new_inttime: u.Quantity):
        """Sets the integration time of the observation.
        Inputs
        - new_inttime float/int or astropy.units.Quantity.
            If no units provided, seconds are assumed.
        """
        if new_inttime <= 0:
            raise ValueError(f"'inttime' must be a positive number, but is {new_inttime}.")

        if not new_inttime.unit.is_equivalent(u.s):
            raise TypeError(f"'inttime' must be a quantity with time units, but is {new_inttime}.")

        self._inttime = new_inttime.to(u.s)

    @property
    def ontarget_fraction(self) -> float:
        """Fraction of the total observing time spent on the target source.
        If scans include more sources than just the target source, then the fraction will be calculated
        from these scans, ignoring this fraction.
        """
        return self._ontarget

    @ontarget_fraction.setter
    @enforce_types
    def ontarget_fraction(self, ontarget: float):
        if not (0.0 < ontarget <= 1.0):
            raise ValueError("'ontarget_fraction' must be a float within (0.0, 1.0].")

        with self._mutex:
            self._ontarget = ontarget
            self._rms = None

    @property
    def ontarget_time(self) -> dict[str, u.Quantity]:
        """Total time spent on the target source during the observation.
        It can return None if the ontarget_fraction and duration have not been set yet, showing a warning.
        """
        if not self.scans:
            return {'DUMMY': self.duration*self.ontarget_fraction}

        to_return = {}
        for scanblock in self.scans.values():
            if len(scanblock.scans) > 1:
                for src, dur in scanblock.fractional_time().items():
                    to_return[src] = dur * self.duration
            else:
                to_return[scanblock.scans[0].source.name] = self.duration*self.ontarget_fraction

        return to_return

    @property
    def bandwidth(self) -> u.Quantity:
        """Returns the total bandwidth of the observation.
        It returns None if the attributes 'polarizations', 'datarate', or 'bitsampling'
        have not been set yet.
        """
        pols = self.polarizations.value % 3 + self.polarizations.value // 3  # Either 1 or 2
        return (self.datarate/(pols*self.bitsampling*2)).to(u.MHz)  # type: ignore

    @property
    def bitsampling(self) -> u.Quantity:
        """Returns the bit sampling at which the data are recorded during the observation.
        """
        return self._bitsampling

    @bitsampling.setter
    @enforce_types
    def bitsampling(self, new_bitsampling: u.Quantity):
        """Sets the bit sampling of the observation.
        Inputs
        - new_bitsampling : int | astropy.units.Quantity (bit-equivalent)
            In bits.
        """
        if not new_bitsampling.unit.is_equivalent(u.bit):
            raise ValueError(f"Unexpected unit for new_bitsampling. Bits spected but {new_bitsampling} received.")

        self._bitsampling = new_bitsampling

    @property
    def stations(self) -> Stations:
        """Returns the network of stations 'Stations' that will participate in this observation
        observing the target source.
        """
        return self._stations

    @stations.setter
    @enforce_types
    def stations(self, new_stations: Stations):
        assert isinstance(new_stations, Stations)
        with self._mutex:
            self._stations = new_stations
            self._elevations = None
            self._altaz = None
            self._is_visible = None
            self._is_always_observable = None
            self._rms = None
            self._baseline_sensitivity = None
            self._uv_baseline = None
            self._uv_array = None
            self._synth_beam = None

    @enforce_types
    def sources_in_block(self, block_name: str, source_type: Optional[SourceType] = None) -> list[Source]:
        """Returns the sources that are included in a specific scan block.
        """
        return self._scans[block_name].sources(source_type)

    @enforce_types
    def sourcenames_in_block(self, block_name: str, source_type: Optional[SourceType] = None) -> list[str]:
        """Returns the sources that are included in a specific scan block.
        """
        return self._scans[block_name].sourcenames(source_type) if self._scans is not None else None

    def elevations(self) -> dict[str, dict[str, coord.angles.Latitude]]:
        """Returns the elevation of the sources for each stations participating in the observation
        for all the given observing times.

        Returns
            elevations : dict
                Dictionary where they keys are: the ScanBlock name, the source name,
                and then the station code names, and the values are the elevations at each time stamp.
        """
        if not self.sources():
            return {}

        with self._mutex:
            if self._elevations is None:
                self._elevations = {src: {ant: altaz.alt for ant, altaz in srcd.items()}
                                    for src, srcd in self.altaz().items()}

        return self._elevations

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
        if not self.sources():
            return {}

        if self._altaz is None:
            self._altaz = self._altaz_threads()

        return self._altaz

    def _altaz_threads(self) -> dict:
        def compute_altaz_for_source(station_source: tuple) -> tuple:
            station, source = station_source[0], station_source[1]
            return station.codename, source.name, station.altaz(self.times, source)

        combined_results: dict[str, dict[str, coord.SkyCoord]] = defaultdict(dict)

        with ThreadPoolExecutor() as pool:
            results = pool.map(compute_altaz_for_source,
                               [(station, source) for station in self.stations
                                for source in self.sources()])
            for station_codename, source_name, altaz in results:
                combined_results[source_name][station_codename] = altaz

        return combined_results

    def _compute_is_visible_for_source(self, ant_src_t: tuple) -> tuple:
        return ant_src_t[0].codename, ant_src_t[1].name, ant_src_t[0].is_observable(ant_src_t[2],
                                                                                    ant_src_t[1])

    @enforce_types
    def is_observable(self, times: Optional[Time] = None) -> dict[str, dict[str, list[bool]]]:
        """Returns whenever the given ScanBlock can be observed by each station for each time
        of the observation. If times are not provided, then it will use the observation times.

        Returns
            is_visible : dict
                Dictionary where they keys are the station code names, and the values will be
                a dictionary with the sources names as keys, and a tuple containing
                a numpy array with the indexes in the `Observation.times`
                array with the times where the target source can be observed by the station.

                In this sense, you can e.g. call obs.times[obs.is_visible[a_station_codename]]
                to get such times.
        """
        if not self.sources():
            return {}

        if self._is_visible is None:
            self._is_visible = {ablockname: {ant.codename: np.full(len(self.times), True, dtype=bool)
                                for ant in self.stations} for ablockname in self.scans}
            with ThreadPoolExecutor() as executor:
                results = list(executor.map(self._compute_is_visible_for_source,
                                            [(station, source, times if times is not None else self.times)
                                             for station in self.stations
                                             for source in self.sources()]))

            for station_codename, source_name, visible in results:
                for ablockname, ablock in self.scans.items():
                    if source_name in ablock.sourcenames():
                        self._is_visible[ablockname][station_codename] &= visible

        return self._is_visible

    def can_be_observed(self) -> dict[str, dict[str, bool]]:
        """Returns whenever the sources can be observed by each station at least during part
        of the observation.

        Returns
            is_observable : dict
                Dictionary where they keys are the scan block names, and the values are
                a dictionary with the station code names as keys, and a boolean containing
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
        if not self.sources():
            return {}

        if self._is_always_visible is None:
            self._is_always_visible = {ablockname: {ant.codename: True
                                       for ant in self.stations} for ablockname in self.scans}

            def compute_is_always_visible_for_source(station_source: tuple) -> tuple:
                station, source = station_source[0], station_source[1]
                return station.codename, source.name, station.is_always_observable(self.times, source)

            with ThreadPoolExecutor() as executor:
                results = list(executor.map(compute_is_always_visible_for_source,
                                            [(station, source) for station in self.stations
                                             for source in self.sources()]))

            for station_codename, source_name, visible in results:
                for ablockname, ablock in self.scans.items():
                    if source_name in ablock.sourcenames():
                        self._is_always_visible[ablockname][station_codename] &= visible

        return self._is_always_visible

    @enforce_types
    def when_is_observable(self, min_stations: int = 3,
                           mandatory_stations: Optional[Literal['all'] | list[str]] = None,
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
            within_time_range = self.times
        else:
            assert len(within_time_range) == 2, \
                "The time range must contain only two times (start and end time)"
            assert within_time_range[1] <= within_time_range[0], \
                "The end time must be larger than the start time"
            within_time_range = within_time_range[0] + \
                np.arange((within_time_range[1] - within_time_range[0]).to(u.day).value, 0.01)*u.day

        # TODO: I can likely optimize this function. And continue from here on
        if mandatory_stations == 'all':
            mandatory_stations = self.stations.station_codenames
            min_stations = len(mandatory_stations)

        result: dict[str, list[Time | u.Quantity]] = {}
        for blockname, station_visibility in self.is_observable(times=within_time_range).items():
            visible_times = [False]*len(within_time_range)
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
        try:
            return ((49500*u.arcsec*u.MHz*u.km)*self.channels /
                    (list(self.longest_baseline().values())[0][1] *
                    self.bandwidth/self.subbands)).to(u.arcsec)
        except IndexError:
            raise SourceNotVisible

    def time_smearing(self) -> u.Quantity:
        """Returns the time smearing expected for the given observation.

        The peak response to a point target source decreases at positions farther away from the
        pointing (correlated) sky position due to the time averaging performed in the data.

        This function returns the angular separation at which the time smearing produces
        a reduction of a 10% in the response of the telescope. The field of view should then
        be limited to this range to avoid significant loses.
        """
        try:
            return ((18560*u.arcsec*u.km*u.s/u.cm) * (self.wavelength /
                    (list(self.longest_baseline().values())[0][1] *
                    self.inttime))).to(u.arcsec)
        except IndexError:
            raise SourceNotVisible

    def datasize(self) -> Optional[u.Quantity]:
        """Returns the expected size for the output FITS IDI files.

        A regular observation with the European VLBI Network is stored in FITS IDI files,
        typically several 2-GB files. This function provides an estimation of the total
        size for these stored files.
        Note that this function does not take into account down times for the different
        stations. The provided value will thus always be un upper-limit for the real, final,
        value.
        """
        temp = len(self.stations)**2 * (self.duration/self.inttime).decompose()
        temp *= self.polarizations.value*self.subbands*self.channels
        return temp * 1.75*u.GB / (131072 * 3600)

    def sun_limiting_epochs(self) -> dict[str, list[Time]]:
        """Returns a dictionary with the sun-constrained epochs for each scan block.
        For each block, it retursn a two-element list with the initial and final time that
        define the time range when the Sun is too close to at least some of the sources
        in the block. The "too close" threshold is defined by the Barray Clark estimates
        from predictions by Ketan Desai of IPM scattering sizes (see SCHED references) at
        the observing band.

        Returns
            dict[str, list[Time]]
                A dictionary with the keys being the name of the scan blocks and the values
                a list with the start and end time for the exclusion time range.
                If the list is empty, implies that there is no risky times where the Sun is too close.
        """
        bad_epochs: dict[str, list[Time]] = {}
        for blockname, block in self.scans.items():
            bad_epochs[blockname] = []
            for src in block.sources():
                bad_epochs[blockname] += src.sun_constraint(freqsetups.min_separation_sun(self.band))

            # Keep only the first and last epochs
            if len(bad_epochs[blockname]) > 0:
                t0, t1 = min(bad_epochs[blockname]), max(bad_epochs[blockname])
                if t0.datetime.month < 3 and t1.datetime.month > 10:
                    mid_year = Time('2025-06-01')
                    t1 = max([t for t in bad_epochs[blockname] if t <= mid_year])
                    t0 = min([t for t in bad_epochs[blockname] if t >= mid_year])

                bad_epochs[blockname] = [t0, t1]

        return bad_epochs

    @enforce_types
    def sun_constraint(self, times: Optional[Time] = None) -> dict[str, Optional[u.Quantity]]:
        """Checks if the Sun is too close to the targets on each block scan, according to the
        Barray Clark estimates from predictions by Ketan Desai of IPM scattering sizes (see
        SCHED references).
        """
        sun_seps: dict[str, Optional[u.Quantity]] = {}
        for blockname, block in self.scans.items():
            min_sep = np.min([s.sun_separation(times=times if times is not None else self._REF_YEAR
                                               if not self.fixed_time else self.times)
                              for s in block.sources()])
            if isinstance(min_sep, float):
                min_sep *= u.deg

            sun_seps[blockname] = min_sep if min_sep <= freqsetups.min_separation_sun(self.band) \
                else None

        return sun_seps

    def thermal_noise(self) -> Union[u.Quantity, dict[str, u.Quantity]]:
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

        # TODO: I think the n_stations = 1 may be a special case, where it requires k=0, not k=1
        sefds = np.array([stat.sefd(self.band).to(u.Jy).value for stat in self.stations])
        bandwidths = np.array([s.datarate.to(u.bit/u.s).value if s.datarate is not None
                              else self.datarate.to(u.bit/u.s).value for s in self.stations]) / \
            (2 * self.bitsampling.to(u.bit).value)  # In Hz
        bandwidth_min = np.minimum.outer(bandwidths, bandwidths)
        if not self.sources():
            if not self.fixed_time:
                dt = self.duration.to(u.s).value * self.ontarget_fraction
            else:
                dt = (self.times[-1]-self.times[0]).to(u.s).value * self.ontarget_fraction

            temp = np.sum((bandwidth_min/np.outer(sefds, sefds)) *
                          np.triu(np.ones_like(bandwidth_min), k=1))
            self._rms = ((1/0.7)/np.sqrt(temp * dt * min(self.polarizations.value, 2)))*u.Jy/u.beam
            return self._rms
        elif not self.fixed_time:
            self._rms = {}
            for sourcename in self.sourcenames:
                temp = np.sum((bandwidth_min/np.outer(sefds, sefds)) *
                              np.triu(np.ones_like(bandwidth_min), k=1))
                self._rms[sourcename] = ((1/0.7)/np.sqrt(temp*self.duration.to(u.s).value *
                                                         min(self.polarizations.value, 2) *
                                                         self.ontarget_fraction))*u.Jy/u.beam

            return self._rms
        else:
            self._rms = {}
            delta_t = (self.times[1] - self.times[0]).to(u.s).value * self.ontarget_fraction
            for sourcename in self.sourcenames:
                scanname = self._scanblock_name_from_source_name(sourcename)
                assert scanname is not None, f"No scan found related to the source {sourcename}"
                visible = np.array([self.is_observable()[scanname][stat.codename] for stat in self.stations])
                integrated_time = np.sum(visible[:, None] & visible[None, :], axis=2) * delta_t
                temp = np.sum((bandwidth_min*integrated_time /
                              np.outer(sefds, sefds)) * np.triu(np.ones_like(bandwidth_min), k=1))
                self._rms[sourcename] = ((1/0.7)/np.sqrt(temp * min(self.polarizations.value, 2)))*u.Jy/u.beam

        return self._rms

    @enforce_types
    def baseline_sensitivity(self, antenna1: Optional[str] = None, antenna2: Optional[str] = None,
                             integration_time: u.Quantity = u.Quantity(1, u.min)) -> u.Quantity:
        """Returns the sensitivity of a given baseline in the observation for
        the expecifiedintegration time.

        If no antennas are defined, then it will return all baseline' sensitivities.
        If only one antenna is specified, then it will return all baseline' sensitivities to that
        particular antenna.

        Inputs
            - antenna1: str | None
                Codename or antenna name of one of the antennas composing the baseline.
            - antenna2: str | None
                Codename or antenna name of the other antennas composing the baseline.
            - integration_time: astropy.units.Quantity  [OPTIONAL]
                Integration time to compute the sensitivity. By default, one minute.
        Returns
            - sensitivity: astropy.units.Quantity
                The sensitivity of the baseline in the observation, in Jy/beam units.
        """
        if self._baseline_sensitivity is None:
            self._baseline_sensitivity = self._baseline_sensitivity_numpy()
            for key in self._baseline_sensitivity:
                self._baseline_sensitivity[key] = self._baseline_sensitivity[key] / \
                    np.sqrt(integration_time.to(u.s).value)
        if antenna1 is None and antenna2 is None:
            return self._baseline_sensitivity
        elif antenna1 is not None and antenna2 is None:
            return [self._baseline_sensitivity[f"{antenna1}-{ant}"]
                    for ant in self.stations.station_codenames] + \
                   [self._baseline_sensitivity[f"{ant}-{antenna1}"]
                    for ant in self.stations.station_codenames if ant != antenna1]
        elif antenna1 is None and antenna2 is not None:
            return [self._baseline_sensitivity[f"{antenna2}-{ant}"]
                    for ant in self.stations.station_codenames] + \
                   [self._baseline_sensitivity[f"{ant}-{antenna2}"]
                    for ant in self.stations.station_codenames if ant != antenna2]

        ant1 = self.stations[antenna1]
        ant2 = self.stations[antenna2]
        if f"{ant1.codename}-{ant2.codename}" in self._baseline_sensitivity:
            return self._baseline_sensitivity[f"{ant1.codename}-{ant2.codename}"]

        return self._baseline_sensitivity[f"{ant2.codename}-{ant1.codename}"]

    # This implementation is ~5 times faster than with threads
    # Only for >20 antennas, it gets ~2 times slower. I think it's worth it.
    def _baseline_sensitivity_numpy(self) -> dict[str, u.Quantity]:
        basel_sens: dict[str, u.Quantity] = {}

        sefds = np.array([stat.sefd(self.band).to(u.Jy).value for stat in self.stations])
        datarates = np.array([s.datarate.to(u.bit/u.s).value if s.datarate is not None
                              else self.datarate.to(u.bit/u.s).value for s in self.stations])
        datarates_min = np.minimum.outer(datarates, datarates)
        sens = (1/0.7)*np.sqrt(np.outer(sefds, sefds)) / \
                                np.sqrt(datarates_min/(min(self.polarizations.value, 2) * \
                                        self.bitsampling.to(u.bit).value))
        for i in range(len(self.stations)):
            for j in range(i, len(self.stations)):
                basel_sens[f"{self.stations[i].codename}-{self.stations[j].codename}"] = \
                    sens[i, j]*u.Jy/u.beam

        return basel_sens

    def _compute_uv_per_source(self, source: Optional[Source] = None) -> dict[str, u.Quantity]:
        nstat = len(self.stations)
        positions = np.array([[xyz.value for xyz in ant.location.to_geocentric()]
                              for ant in self.stations])

        # Compute all unique baselines (i < j)
        idx_i, idx_j = np.triu_indices(nstat, k=1)
        bl_xyz = positions[idx_i] - positions[idx_j]  # shape (nbl, 3)
        bl_names = [f"{self.stations[int(i)].codename}-{self.stations[int(j)].codename}"
                    for i, j in zip(idx_i, idx_j)]

        gstimes = self.gstimes
        # Prepare the rotation matrix for each time
        if source is None:
            hourangle: coord.Angle = gstimes
            hourangle = gstimes
            dec = 45 * u.deg
            # print:"Target source is not set, thus we assume a source at +/- 45º declination"
            #       " to estimate the (u, v) values.'")
        else:
            hourangle = (gstimes - source.ra.to(u.hourangle)).wrap_at(24 * u.hourangle)
            dec = source.dec

        ha_rad = hourangle.to(u.rad).value  # shape (ntimes,)
        dec_rad = dec.to(u.rad).value

        # Rotation matrix for each time (ntimes, 3, 3)
        # m = np.zeros((len(ha_rad), 3, 3))
        # m[:, 0, 0] = np.sin(ha_rad)
        # m[:, 0, 1] = np.cos(ha_rad)
        # # m[:, 0, 2] = 0
        # m[:, 1, 0] = -np.sin(dec_rad) * np.cos(ha_rad)
        # m[:, 1, 1] = np.sin(dec_rad) * np.sin(ha_rad)
        # m[:, 1, 2] = np.cos(dec_rad)
        # # m[:, 2, :] = 0  # If you need w, fill in as needed

        # m = np.array([[np.sin(ha_rad), np.cos(ha_rad), np.zeros(len(ha_rad))],
        #               [-np.sin(dec_rad)*np.cos(ha_rad),
        #                np.sin(dec_rad)*np.sin(ha_rad),
        #                np.cos(dec_rad)*np.ones(len(ha_rad))], np.zeros((3, len(ha_rad)))])
        m = np.array([[[np.sin(h), np.cos(h), 0],
                       [-np.sin(dec_rad)*np.cos(h), np.sin(dec_rad)*np.sin(h), np.cos(dec_rad)],
                       [0, 0, 0]] for h in ha_rad])

        # Compute uvw for all baselines and times: (ntimes, nbl, 3)
        # Use einsum for batch matrix multiplication
        bl_uv = (np.einsum('tij,bj->tbi', m, bl_xyz)*u.m/self.wavelength).decompose()
        # shape (ntimes, nbl, 3)
        bl_uv_up: dict[str, u.Quantity] = {}
        if source is None:
            for i, bl_name in enumerate(bl_names):
                bl_uv_up[bl_name] = bl_uv[:, i, :2]  # Only u,v
            return bl_uv_up

        # Source is provided: filter by visibility
        ants_up = self.is_observable()
        blockname = [ablock for ablock in self.scans if source.name in self.scans[ablock].sourcenames()][0]

        for bl_idx, bl_name in enumerate(bl_names):
            ant1, ant2 = bl_name.split('-')
            # Find times where both antennas are up
            # up_times = np.intersect1d(ants_up[blockname][ant1], ants_up[blockname][ant2],
            #                           assume_unique=True)
            # NOTE: type(ants_up)=<class 'dict'>, type(ants_up[blockname])=<class 'dict'>,
            # type(ants_up[blockname][ant1])=<class 'numpy.ndarray'>
            up_times = ants_up[blockname][ant1] & ants_up[blockname][ant2]  # type: ignore
            if len(up_times) > 0:
                # up_times are indices into the time array
                bl_uv_up[bl_name] = bl_uv[up_times, bl_idx, :2]

        if not bl_uv_up:
            raise SourceNotVisible

        return bl_uv_up

    def _compute_uv_per_source_legacy(self, source: Optional[Source] = None) -> dict[str, u.Quantity]:
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

        gstimes = self.gstimes
        if source is None:
            hourangle: coord.Angle = gstimes
            print("WARNING: 'target' is not set, thus we assume a source at +/- 45º declination"
                  " to estimate the (u, v) values.'")
            m = np.array([[np.sin(hourangle), np.cos(hourangle), np.zeros(len(hourangle))],
                          [-np.sin(45*u.deg)*np.cos(hourangle),
                           np.sin(45*u.deg)*np.sin(hourangle),
                           np.cos(45*u.deg)*np.ones(len(hourangle))]])
        else:
            hourangle = (gstimes - source.ra.to(u.hourangle)).value % 24*u.hourangle
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
        if not self.sources():
            return {'DUMMY': self._compute_uv_per_source()}

        with self._mutex_uv:
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

        def resolution(bl: float) -> u.Quantity:
            return ((2.063e8*u.mas)/bl).to(u.mas)

        uvvis = self.get_uv_values()
        for src, uv in uvvis.items():
            # Transform the uv points into r,theta (polar) points
            uvvis_polar = np.empty_like(uv)
            uvvis_polar[:, 0] = np.sqrt((uv**2).sum(axis=1))  # radius
            uvvis_polar[:, 1] = np.arctan2(uv[:, 1], uv[:, 0])  # theta
            # Defines the BMAJ and PA
            bl_bmaj = np.max(uvvis_polar[:, 0])
            bl_bmaj_theta = uvvis_polar[:, 1][np.where(uvvis_polar[:, 0] == bl_bmaj)][0]
            # Gets the BMIN and an orthogonal projection
            bl_bmin_theta = (bl_bmaj_theta + np.pi/2) % (2*np.pi)
            bl_bmin = np.max(np.abs(uv.dot(np.array([np.cos(bl_bmin_theta),
                                                     np.sin(bl_bmin_theta)]))))

            self._synth_beam[src] = {'bmaj': resolution(bl_bmin), 'bmin': resolution(bl_bmaj),
                                     'pa': (bl_bmaj_theta*u.rad).to(u.deg)}

        return self._synth_beam

    @enforce_types
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
        raise NotImplementedError
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
            uvimg[oversampling*pixsize//2 + np.trunc(uvdata[:, 0]*uvscaling).astype(int),
                  oversampling*pixsize//2 + np.trunc(uvdata[:, 1]*uvscaling).astype(int)] += 1
        else:
            uvimg[oversampling*pixsize//2 + np.trunc(uvdata[:, 0]*uvscaling).astype(int),
                  oversampling*pixsize//2 + np.trunc(uvdata[:, 1]*uvscaling).astype(int)] = 1

        # Recovering the requested size for the image
        # uvimg = uvimg[int(pixsize*0.5):int(pixsize*1.5), int(pixsize*0.5):int(pixsize*1.5)]
        dirty_beam = np.real(np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(uvimg/np.max(uvimg)))))
        imgsize = (uvscaling*u.rad/2).to(u.mas)  # angular equivalent size of the resulting image
        return dirty_beam[int(pixsize*(oversampling-1)/2):int(pixsize*(oversampling+1)/2),
                          int(pixsize*(oversampling-1)/2):int(pixsize*(oversampling+1)/2)].T / \
            np.max(dirty_beam), np.linspace(-imgsize, imgsize, pixsize)

    @enforce_types
    def print_obs_times(self, date_format: str = '%d %B %Y'):
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

    # def get_fig_dirty_map(self):
    #     raise NotImplementedError
    #     # Right now I only leave the natural weighting map (the uniform does not always correspond to the true one)
    #     dirty_map_nat, laxis = self.get_dirtymap(pixsize=1024, robust='natural', oversampling=4)
    #     # TODO: uncomment these two lines and move them outside observation. Flexibility
    #     # fig1 = px.imshow(img=dirty_map_nat, x=laxis, y=laxis[::-1], labels={'x': 'RA (mas)', 'y': 'Dec (mas)'}, \
    #     #         aspect='equal')
    #     # fig = make_subplots(rows=1, cols=1, subplot_titles=('Natural weighting',),
    #     #                     shared_xaxes=True, shared_yaxes=True)
    #     fig.add_trace(fig1.data[0], row=1, col=1)
    #     mapsize = 30*self.synthesized_beam()['bmaj'].to(u.mas).value
    #     fig.update_layout(coloraxis={'showscale': False, 'colorscale': 'Inferno'}, showlegend=False,
    #                       xaxis={'autorange': False, 'range': [mapsize, -mapsize]},
    #                       yaxis={'autorange': False, 'range': [-mapsize, mapsize]}, autosize=False)
    #     fig.update_xaxes(title_text="RA (mas)", constrain="domain")
    #     fig.update_yaxes(title_text="Dec (mas)", scaleanchor="x", scaleratio=1)
    #     # dirty_map_nat, laxis = obs.get_dirtymap(pixsize=1024, robust='natural', oversampling=4)
    #     # dirty_map_uni, laxis = obs.get_dirtymap(pixsize=1024, robust='uniform', oversampling=4)
    #     # fig1 = px.imshow(img=dirty_map_nat, x=laxis, y=laxis[::-1], labels={'x': 'RA (mas)', 'y': 'Dec (mas)'}, \
    #     #         aspect='equal')
    #     # fig2 = px.imshow(img=dirty_map_uni, x=laxis, y=laxis[::-1], labels={'x': 'RA (mas)', 'y': 'Dec (mas)'}, \
    #     #         aspect='equal')
    #     # fig = make_subplots(rows=1, cols=2, subplot_titles=('Natural weighting', 'Uniform weighting'),
    #     #                     shared_xaxes=True, shared_yaxes=True)
    #     # fig.add_trace(fig1.data[0], row=1, col=1)
    #     # fig.add_trace(fig2.data[0], row=1, col=2)
    #     # mapsize = 30*obs.synthesized_beam()['bmaj'].to(u.mas).value
    #     # fig.update_layout(coloraxis={'showscale': False, 'colorscale': 'Inferno'}, showlegend=False,
    #     #                   xaxis={'autorange': False, 'range': [mapsize, -mapsize]},
    #     #                   # This xaxis2 represents the xaxis for fig2.
    #     #                   xaxis2={'autorange': False, 'range': [mapsize, -mapsize]},
    #     #                   yaxis={'autorange': False, 'range': [-mapsize, mapsize]}, autosize=False)
    #     # fig.update_xaxes(title_text="RA (mas)", constrain="domain")
    #     # fig.update_yaxes(title_text="Dec (mas)", row=1, col=1, scaleanchor="x", scaleratio=1)
    #
    #     return fig

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
