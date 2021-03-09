# -*- coding: utf-8 -*-
# Licensed under GPLv3+ - see LICENSE
from __future__ import annotations
import configparser
from importlib import resources
import numpy as np
from astropy import units as u
from astropy import coordinates as coord
from astropy.io import ascii
from astropy.time import Time
from astroplan import Observer, FixedTarget

"""Module that defines the `Station` and `Stations` objects, which represent a station (antenna)
or a network composed of antennas.
"""

__all__ = ['Station', 'SelectedStation', 'Stations']



class Station(object):
    """Defines an astronomical station (antenna).
    A station is defined by some names and codenames, coordinates, and its sensitivity for the
    radio bands (wavelengths) it can observe.
    Apart of the metadata related to the station, it allows to compute the altitude/azimuth, elevation,
    or simply when a source is visible from the station for a given time range.
    """
    def __init__(self, name: str, codename: str, network: str, location: coord.EarthLocation,
                 freqs_sefds: dict, min_elevation=20*u.deg, fullname: str = None,
                 all_networks: str = None, country: str = '', diameter: str = '', real_time: bool = False):
        """Initializes a station.

        Inputs
        - name : str
            Name of the observer (the station that is going to observe).
            If it contains undercores (_), they will be converted to blank spaces.
        - codename : str
            A short code (accronym) for the name of the station. It is meant to follow the standard approach
            from the EVN: an (often) two-letter code unique for each station.
        - network : str
            Name of the network to which the station belongs (e.g. EVN).
        - location : astropy.coordinates.EarthLocation
            Position of the station on Earth in (x,y,z) gecentric coordinates.
        - freqs_sefds : dict
            Dictionary with all frequencies the station can observe as keys of the dictionary, and the
            values representing the system equivalent flux density (SEFD; in Jansky units)
            at each frequency.
            Although the key format is in principle free, we recommend to use the syntax 'XXcm' (str type).
            This will be then consistent with the default station catalog.
        - min_elevation : Quantity [OPTIONAL]
            Minimum elevation that the station can reach to observe a source. If no units (astropy.units)
            provided, degrees are assumed. By default it 20 degrees. It does not support an azimuth-dependent
            elevation limits.
        - fullname : str [OPTIONAL]
            Full name of the station. If not given, same as `name` is assumed.
            It can be used to expand the full name if an abbreviation is typically used for the name.
            For example, name: VLA, fullname: Karl G. Jansky Very Large Array.
        - all_networks : str [OPTIONAL]
            Networks where the station can participate (free style string).
        - country : str [OPTIONAL]
            Country where the station is located.
        - diameter : str [OPTIONAL]
            Diameter of the station (free format string). We recommend a syntax of e.g. '30 m' for normal
            single-dish antennas, and in case of interferometers it can have a form like '25 x 20 m',
            meaning that the station is composed of 25 antennas of 20 m each.
        - real_time : bool [OPTIONAL]
            If the station can participate in real-time observations (e.g. e-EVN), False by default.
        """
        self.observer = Observer(name=name.replace('_', ' '), location=location)
        self._codename = codename
        self._network = network
        self._freqs_sefds = freqs_sefds
        if (type(min_elevation) is float) or (type(min_elevation) is int):
            self._min_elev = min_elevation*u.deg
        else:
            self._min_elev = min_elevation

        if fullname is None:
            self._fullname = name
        else:
            self._fullname = fullname

        if all_networks is None:
            self._all_networks = network
        else:
            self._all_networks = all_networks

        self._country = country
        self._diameter = diameter
        self._real_time = real_time


    @property
    def name(self) -> str:
        """Name of the station.
        """
        return self.observer.name


    @property
    def codename(self) -> str:
        """Codename of the station (typically a two-letter accronym).
        """
        return self._codename


    @property
    def fullname(self) -> str:
        """Full name of the station. If not specified, it can be the same as 'name'.
        """
        return self._fullname


    @property
    def network(self) -> str:
        """Name of the network to which the station belongs.
        """
        return self._network


    @property
    def all_networks(self) -> str:
        """Name of all networks to which the station belongs.
        If not specified it can be the same as 'network'.
        """
        return self._all_networks


    @property
    def country(self) -> str:
        """Country where this station is located.
        It can be an empty string if not specified.
        """
        return self._country


    @property
    def diameter(self) -> str:
        """String representing the diameter of the station, and/or how many antennas compose
        the station in case of connected-interferometers.
        """
        return self._diameter


    @property
    def real_time(self) -> bool:
        """If the station can participate in real-time observations (e.g. e-EVN).
        """
        return self._real_time


    @property
    def location(self) -> coord.EarthLocation:
        """Location of the station as an astropy.coordinates.EarthLocation object.
        """
        return self.observer.location


    @property
    def bands(self):
        """Observing bands the station can observe.
        Returns a dict_keys object with all bands in a string format as introduced in the freqs_sefd
        attribute when the Station was created.
        """
        return self._freqs_sefds.keys()


    @property
    def sefds(self) -> dict:
        """Returns a dictionary with the system equivalent flux density (SEFDs) for each
        of the frequencies the station can observe (given as keys).
        """
        return self._freqs_sefds


    @property
    def min_elevation(self) -> u.Quantity:
        """Minimum elevation the station can observe a source.
        Returns an astropy.units.Quantity (i.e. number with units).
        """
        return self._min_elev


    def elevation(self, obs_times: Time, target: FixedTarget) -> coord.angles.Latitude:
        """Returns the elevation of the target source as seen by the Station during obs_times.

        Inputs
        - obs_times : astropy.time.Time
            Time to compute the elevation of the source
            (either a single timestamp or an array of times).
        - target : astroplan.FixedTarget
             Target to observe.

        Output
        - elevations : astropy.coordinates.angles.Latitute
            Elevation of the source at the given obs_times.
        """
        return self.observer.altaz(obs_times, target).alt


    def altaz(self, obs_times: Time, target: FixedTarget) -> coord.sky_coordinate.SkyCoord:
        """Returns the altaz coordinates of the target source for the given observing times.

        Inputs
        - obs_times : astropy.time.Time
            Time to compute the elevation of the source
            (either a single timestamp or an array of times).
        - target : astroplan.FixedTarget
             Target coordinates to observe.

        Output
        - altaz : astropy.coordinates.sky_coordinate.SkyCoord
            Altitude and azimuth of the source for each given time.
        """
        return self.observer.altaz(obs_times, target)


    def is_visible(self, obs_times: Time, target: FixedTarget) -> tuple:
        """Returns when the target source is visible for this station at the given times.

        Inputs
        - obs_times : astropy.time.Time
            Time to compute the elevation of the source
            (either a single timestamp or an array of times).
        - target : astroplan.FixedTarget
             Target coordinates to observe.

        Output
        - visible : tuple
            Tuple containing the indexes of obs_times when the target source is visible
            from the station. Therefore obs_times[visible] would return only those times.
        """
        elevations = self.elevation(obs_times, target)
        return np.where(elevations >= self.min_elevation)


    def has_band(self, band: str) -> bool:
        """Returns if the Station can observed the given band `the_band`.

        Inputs
        - band : str
            A string representing an observing band, following the same syntax as used
            when the station was initialized and the bands where defined in the keys of
            the freqs_sefds attribute.

        Output
        - bool whenever the station has the given observing band.
        """
        return band in self.bands


    def sefd(self, band: str) -> float:
        """Returns the system equivalent flux density (SEFD) of the Station at the given band,
        in Jansky (Jy) units.

        Input
        - band : str
            A string representing an observing band, following the same syntax as used
            when the station was initialized and the bands where defined in the keys of
            the freqs_sefds attribute.

        Output
        - SEFD : float
            The SEFD at the given band, in Jy units.

        Exception
        - It may raise KeyError if the given band is not available for this station.

        """
        return self._freqs_sefds[band]


    def __str__(self):
        return f"<{self.codename}>"


    def __repr__(self):
        return f"<Station: {self.codename}>"



class SelectedStation(Station):
    """Extends the Station class with an additional attribute: `selected` (bool).
    This allows the sub-seting of lists of stations by considering all of them to observe
    or just disabling some of them.
    """
    def __init__(self, name: str, codename: str, network: str, location: coord.EarthLocation,
                 freqs_sefds: dict, min_elevation=20*u.deg, fullname: str = None,
                 all_networks: str = None, country: str = '', diameter: str = '',
                 real_time: bool = False, selected: bool = True):
        """Initializes a SelectedStation.
        This class extends Station by adding one additional attribute: `selected` (bool).

        Inputs
        - name : str
            Name of the observer (the station that is going to observe).
            If it contains undercores (_), they will be converted to blank spaces.
        - codename : str
            A short code (accronym) for the name of the station. It is meant to follow the standard approach
            from the EVN: an (often) two-letter code unique for each station.
        - network : str
            Name of the network to which the station belongs (e.g. EVN).
        - location : astropy.coordinates.EarthLocation
            Position of the station on Earth in (x,y,z) gecentric coordinates.
        - freqs_sefds : dict
            Dictionary with all frequencies the station can observe as keys of the dictionary, and the
            values representing the system equivalent flux density (SEFD; in Jansky units)
            at each frequency.
            Although the key format is in principle free, we recommend to use the syntax 'XXcm' (str type).
            This will be then consistent with the default station catalog.
        - min_elevation : Quantity [OPTIONAL]
            Minimum elevation that the station can reach to observe a source. If no units (astropy.units)
            provided, degrees are assumed. By default it 20 degrees. It does not support an azimuth-dependent
            elevation limits.
        - fullname : str [OPTIONAL]
            Full name of the station. If not given, same as `name` is assumed.
            It can be used to expand the full name if an abbreviation is typically used for the name.
            For example, name: VLA, fullname: Karl G. Jansky Very Large Array.
        - all_networks : str [OPTIONAL]
            Networks where the station can participate (free style string).
        - country : str [OPTIONAL]
            Country where the station is located.
        - diameter : str [OPTIONAL]
            Diameter of the station (free format string). We recommend a syntax of e.g. '30 m' for normal
            single-dish antennas, and in case of interferometers it can have a form like '25 x 20 m',
            meaning that the station is composed of 25 antennas of 20 m each.
        - real_time : bool [OPTIONAL]
            If the station can participate in real-time observations (e.g. e-EVN), False by default.
        - selected : bool [OPTIONAL]
            If the station is selected to participate in a given observation or not. True by default.
        """

        self._selected = selected
        super().__init__(name, codename, network, location, freqs_sefds,
                         min_elevation, fullname, all_networks, country, diameter, real_time)


    @property
    def selected(self) -> bool:
        """If the station is selected to participate in a given observation or not.
        If False, it will not observe.
        """
        return self._selected


    @selected.setter
    def selected(self, isselected: bool):
        assert isinstance(isselected, bool)
        self._selected = isselected



class Stations(object):
    """Defines a network (collection) of stations (`Station` objects) that
    can participate in an observation together.
    """
    def __init__(self, name: str, stations: list):
        """Initializes a Stations.

        Inputs
        - name : str
            Name associated to the network of stations.
        - stations : list of Station-type elements
            List with all stations belonging to the given network.
        """
        self._name = name
        self._stations = {}
        for a_station in stations:
            if a_station.codename not in self._stations.keys():
                self._stations[a_station.codename] = a_station

        self._codenames = tuple(self._stations.keys())


    @property
    def name(self) -> str:
        """Name of the network of stations.
        """
        return self._name

    @name.setter
    def name(self, new_name: str):
        """Assigns a new name to the network of stations
        """
        self._name = new_name


    @property
    def stations(self) -> list:
        """Returns a list containing all stations in the network.
        """
        return list(self._stations.values())


    @property
    def number_of_stations(self) -> int:
        """Returns the total number of stations in the network.
        """
        return len(self.stations)


    @property
    def codenames(self) -> tuple:
        """Returns a tuple with the `codenames` from all the stations in the network.
        """
        return self._codenames


    @property
    def observing_bands(self) -> set:
        """Returns a set with all `bands` that the stations in the network can observe,
        or at least a subset of stations.
        """
        bands = set()
        for a_station in self._stations.values():
                station_bands = set(a_station.bands)
                bands.update(station_bands)

        return bands



    def add(self, a_station: Station):
        """Adds a new station to the network.
        If a station with the same codename is already present, it will do nothing.

        Inputs
        - a_station : Station
            Station to be added to the network.
        """
        assert isinstance(a_station, Station)
        if a_station.codename in self.codenames:
            print(f"WARNING: {a_station.codename} already in {self.name}. Ignoring addition.")
        else:
            self._stations[a_station.codename] = a_station
            self._codenames = tuple(self._stations.keys())


    def __str__(self):
        return f"<{self.name}: <{', '.join(self.codenames)}>>"


    def __len__(self):
        return self._stations.__len__()


    def __getitem__(self, key):
        if isinstance(key, int):
            return self._stations[self.codenames[key]]
        else:
            return self._stations[key]


    def __setitem__(self, key, value):
        self._stations[key] = value
        self._codenames = tuple(self._stations.keys())


    def __delitem__(self, key):
        if isinstance(key, int):
            self._stations.__delitem__(self.keys[key])
        else:
            self._stations.__delitem__(key)

        self._codenames = tuple(self._stations.keys())


    def __iter__(self):
        return iter(self._stations.values())


    def __contains__(self, item):
        return self._stations.__contains__(item)


    @staticmethod
    def get_stations_from_configfile(filename: str = None, codenames: list = None,
                                               name: str = 'network') -> Stations:
        """Creates a Stations object (i.e. a network of stations) by reading the station
        information from an input file. Optionally, it allows to select only a subset of
        all stations in the file.

        Inputs
        - filename : str
            Path to the text file containing the information from all stations.
            If not provided, it reads the default station catalog file located in
                            data/stations_catalog.inp
            Any other file should contain the same format (standard Python input config files),
            with the following fields per station (whose name would be provided as the name of
            section).
            - station - full name of the station.
            - code : codename assigned to the station. It must be unique (typically two letters).
            - network - main network to which it belongs to.
            - possible_networks - all networks the station can participate in (including 'network')
            - country - country where the station is located.
            - diameter - free format string with the diameter of the station
                        (optional more information in case of interferometers).
            - position = x, y, z (in meters). Geocentric position of the station.
            - min_elevation (in degrees) - minimum elevation the station can observe.
            - real_time = yes/no - if the station can participate in real-time observations (e.g. e-EVN).
            - SEFD_**  - SEFD (in Jy units) of the station at the **cm band. If a given band is not present,
                        it is assumed that the station cannot observe it.
                        For example SEFD_21 = 500 means that the SEFD at 21cm is 500 Jy.
            - Any other attribute is accepted, but ignored in this code. That would easily allow future
              extensions of the code.
        - codenames : list
            If you only want to select a subset of all stations available in the input file,
            here you can pass a list with the codenames of the stations that should be imported.
        - name : str
            Name to assign to the network of stations that will be created.

        Returns
        - network : Stations
            Returns a Stations object containing the selected stations.
        """
        config = configparser.ConfigParser()
        if filename is None:
            with resources.path("data", "stations_catalog.inp") as stations_catalog_path:
                config.read(stations_catalog_path)
        else:
            # With this approach it raises a FileNotFound exception.
            # Otherwise config will run smoothly and provide an empty list.
            config.read(open(filename, 'r'))

        networks = Stations(name, [])
        for stationname in config.sections():
            if (codenames is None) or (config[stationname]['code'] in codenames):
                temp = [float(i.strip()) for i in config[stationname]['position'].split(',')]
                a_loc = coord.EarthLocation(temp[0]*u.m, temp[1]*u.m, temp[2]*u.m)
                # Getting the SEFD values for the bands
                min_elev = float(config[stationname]['min_elevation'])*u.deg
                does_real_time = True if config[stationname]['real_time']=='yes' else False
                sefds = {}
                for akey in config[stationname].keys():
                    if 'SEFD_' in akey.upper():
                        sefds[f"{akey.upper().replace('SEFD_', '').strip()}cm"] = \
                                            float(config[stationname][akey])

                new_station = SelectedStation(stationname, config[stationname]['code'],
                        config[stationname]['network'], a_loc, sefds, min_elev,
                        config[stationname]['station'], config[stationname]['possible_networks'],
                        config[stationname]['country'], config[stationname]['diameter'], does_real_time)
                networks.add(new_station)

        return networks


    def stations_with_band(self, band: str, output_network_name: str = None) -> Stations:
        """Given the current network, it creates a sub-network including only the stations
        that can observe at the given band.

        Inputs
        - band : str
            The observing band that will be selected in the stations.
        - output_network_name : str [OPTIONAL]
            The name assigned to the new network.
            By default, it will be the same as the original network followed by @{band}.

        Returns
        - subnetwork : Stations
            A subset of the original Stations containing only the stations that can observe
            at the given band.
        """
        if output_network_name is None:
            output_network_name = f"{self.name}@{band}"

        subnetwork = Stations(output_network_name, [])
        for station in self.stations:
            if band in station.bands:
                subnetwork.add(station)

        return subnetwork


    def select_stations(self, codenames: list, name: str = None) -> Stations:
        """Returns a new Stations object which will only contain the stations
        defined by the given list of codenames. It will thus be a subset of the current
        network.

        Input
        - codenames : list
            List with the codenames of the stations that should be present in the new
            network.
        - name : str [OPTIONAL]
            Name assigned to the new network. If not provided, it will be the original
            name with the 'sub' preffix.
        Returns
        - subnetwork : Stations
            A new Stations object containing only the defined stations.
        Exceptions
        - It may raise KeyError if one of the given codenames are not present
          among the current stations.
        """
        subnetwork = Stations(name if name is not None else f"sub{self.name}", [])
        for codename in codenames:
            subnetwork.add(self.stations[codename])

        return subnetwork


    @staticmethod
    def get_network_names_from_configfile(filename: str) -> dict:
        """Reads a config file containing the different VLBI networks defined as a config parser file.
        Returns a dictionary with the nickname of the VLBI network as keys,
        Inputs
        - filename : str
            Path to the text file containing the information from all stations.
            If not provided, it reads the default station catalog file located in
                            data/stations_catalog.inp
            Any other file should contain the same format (standard Python input config files),
            with the following fields per station (whose name would be provided as the name of
            section).
            - station - full name of the station.
            - code : codename assigned to the station. It must be unique (typically two letters).
            - network - main network to which it belongs to.
            - possible_networks - all networks the station can participate in (including 'network')
            - country - country where the station is located.
            - diameter - free format string with the diameter of the station
                        (optional more information in case of interferometers).
            - position = x, y, z (in meters). Geocentric position of the station.
            - min_elevation (in degrees) - minimum elevation the station can observe.
            - real_time = yes/no - if the station can participate in real-time observations (e.g. e-EVN).
            - SEFD_**  - SEFD (in Jy units) of the station at the **cm band. If a given band is not present,
                        it is assumed that the station cannot observe it.
                        For example SEFD_21 = 500 means that the SEFD at 21cm is 500 Jy.
            - Any other attribute is accepted, but ignored in this code. That would easily allow future
              extensions of the code.
        - codenames : list
            If you only want to select a subset of all stations available in the input file,
            here you can pass a list with the codenames of the stations that should be imported.
        - name : str
            Name to assign to the network of stations that will be created.

        Returns
        - network : Stations
            Returns a Stations object containing the selected stations.
        """
        config = configparser.ConfigParser()
        if filename is None:
            with resources.path("data", "stations_catalog.inp") as stations_catalog_path:
                config.read(stations_catalog_path)
        else:
            # With this approach it raises a FileNotFound exception.
            # Otherwise config will run smoothly and provide an empty list.
            config.read(open(filename, 'r'))

        networks = Stations(name, [])
        for stationname in config.sections():
            if (codenames is None) or (config[stationname]['code'] in codenames):
                temp = [float(i.strip()) for i in config[stationname]['position'].split(',')]
                a_loc = coord.EarthLocation(temp[0]*u.m, temp[1]*u.m, temp[2]*u.m)
                # Getting the SEFD values for the bands
                min_elev = float(config[stationname]['min_elevation'])*u.deg
                does_real_time = True if config[stationname]['real_time']=='yes' else False
                sefds = {}
                for akey in config[stationname].keys():
                    if 'SEFD_' in akey.upper():
                        sefds[f"{akey.upper().replace('SEFD_', '').strip()}cm"] = \
                                            float(config[stationname][akey])

                new_station = SelectedStation(stationname, config[stationname]['code'],
                        config[stationname]['network'], a_loc, sefds, min_elev,
                        config[stationname]['station'], config[stationname]['possible_networks'],
                        config[stationname]['country'], config[stationname]['diameter'], does_real_time)
                networks.add(new_station)

        return networks




