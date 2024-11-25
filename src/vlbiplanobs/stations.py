# -*- coding: utf-8 -*-
# Licensed under GPLv3+ - see LICENSE
from __future__ import annotations
from collections import abc
from typing import Optional, Union, Iterable, Sequence, Self
import configparser
from importlib import resources
import numpy as np
from rich import print as rprint
from astropy import units as u
from astropy import coordinates as coord
# from astropy.io import ascii
from astropy.time import Time
from astroplan import Observer, FixedTarget, is_observable, is_always_observable
from dataclasses import dataclass
from enum import Enum, auto
from . import constraints
from . import freqsetups

"""Module that defines the `Station` and `Network` objects, which represent a station (antenna)
or a network composed of antennas.
"""

__all__: list[str] = ['Station', 'SelectedStation', 'Network', 'Mount', 'MountType']

class MountType(Enum):
    ALTAZ = auto()
    EQUAT = auto()


@dataclass
class Axis:
    limits: tuple[u.Quantity, u.Quantity]
    speed: u.Quantity
    acceleration: u.Quantity = u.Quantity(0.0, u.deg/u.s/u.s)


@dataclass
class Mount:
    """Defines a telescope mount, including the mount type (ALTAZ, EQUAT), the limits of the telescope
    to point in the sky, and the slewing velocity (rate) and acceleration (acc).
    """
    mount_type: MountType
    ax1: Axis
    ax2: Axis


class Station(object):
    """Defines an astronomical station (antenna).
    A station is defined by some names and codenames, coordinates, and its sensitivity for the
    radio bands (wavelengths) it can observe.
    Apart of the metadata related to the station, it allows to compute the altitude/azimuth, elevation,
    or simply when a source is visible from the station for a given time range.
    """
    def __init__(self, name: str, codename: str, network: str, location: coord.EarthLocation,
                 freqs_sefds: dict[str, u.Quantity],
                 fullname: Optional[str] = None, all_networks: Optional[list[str]] = None, country: str = '',
                 diameter: str = '', real_time: bool = False, mount: Optional[Mount] = None,
                 max_datarate: Optional[u.Quantity | dict[str, u.Quantity]] = None) -> None:
        """Initializes a station.

        Inputs
        - name : str
            Name of the station.
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
        - fullname : str [OPTIONAL]
            Full name of the station. If not given, same as `name` is assumed.
            It can be used to expand the full name if an abbreviation is typically used for the name.
            For example, name: VLA, fullname: Karl G. Jansky Very Large Array.
        - all_networks : list[str] [OPTIONAL]
            All networks where the station can participate.
            By default it will always include 'network' too.
        - country : str [OPTIONAL]
            Country where the station is located.
        - diameter : str [OPTIONAL]
            Diameter of the station (free format string). We recommend a syntax of e.g. '30 m' for normal
            single-dish antennas, and in case of interferometers it can have a form like '25 x 20 m',
            meaning that the station is composed of 25 antennas of 20 m each.
        - real_time : bool [OPTIONAL]
            If the station can participate in real-time observations (e.g. e-EVN), False by default.
        - mount : Mount [OPTIONAL]
            The mount of the station, including the type of mount and the slewing limits, speed, and acceleration.
            If not provided, it will assume an ALTAZ mount with no pointing limits and very high slewing speed.
        - max_datarate : u.Quantity or dict[str, u.Quantity]  [OPTIONAL]
            Specifies the maximum data rate that the station can record. If not provided, it will use the
            data rate assumed in the observation. It can be either a astropy.unit.Quantity value
            (equivalent to Mb/s), or a dictionary in the case that the antenna can show different data rates
            when it is participating in different networks. In that case, the abbreviation of the network will
            be the key of the dictionary, with the quantity as value. For example, Darnhall can record at 4 Gbps
            when participating within e-MERLIN observations, but their data rate is limited to 512 Mbps when
            participating within the EVN.
        """
        # Some sanity checks
        for a_var, a_var_name in zip((name, codename, network, country, diameter), \
                                     ('name', 'codename', 'network', 'country', 'diameter')):
            assert isinstance(a_var, str), f"'{a_var_name}' must be a str."

        assert isinstance(fullname, str) or fullname is None
        assert isinstance(real_time, bool), "'real_time' must be a bool."
        self.observer: Observer = Observer(name=name.replace('_', ' '), location=location)
        self._codename: str = codename
        self._network: str = network
        self._freqs_sefds: dict[str, float] = {f if 'cm' in f else f"{f}cm": v for f,v in freqs_sefds.items()}
        self._fullname: str = name if fullname is None else fullname
        assert isinstance(all_networks, Sequence) or all_networks is None
        self._all_networks: list[str] = [network] if all_networks is None else list(all_networks)
        if network not in self._all_networks:
            self._all_networks.insert(0, network)

        self._country: str = country
        self._diameter: str = diameter
        self._real_time: bool = real_time
        if mount is None:
            self._mount: Mount = Mount(MountType.ALTAZ,
                                       Axis((-10*u.deg, 370*u.deg), 300*u.deg/u.s, 0.0*u.deg/u.s/u.s),
                                       Axis((10*u.deg, 100*u.deg), 300*u.deg/u.s, 0.0*u.deg/u.s/u.s))
        else:
            assert isinstance(mount, Mount)
            self._mount = mount

        self._max_datarate: u.Quantity | dict | None = max_datarate
        if self.mount.mount_type.ALTAZ:
            self._constraints = [constraints.AzimuthConstraint(min=self.mount.ax1.limits[0],
                                                                 max=self.mount.ax1.limits[1]),
                                 constraints.ElevationConstraint(min=self.mount.ax2.limits[0],
                                                                   max=self.mount.ax2.limits[1])]
        else:
            self._constraints = [constraints.HourAngleConstraint(min=self.mount.ax1.limits[0],
                                                                 max=self.mount.ax1.limits[1]),
                                 constraints.DeclinationConstraint(min=self.mount.ax2.limits[0],
                                                                   max=self.mount.ax2.limits[1])]


    @property
    def name(self) -> str:
        """Name of the station.
        """
        assert self.observer.name is not None
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
    def all_networks(self) -> list[str]:
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
    def constraints(self) -> list[constraints.Constraint]:
        """Returns the observing constraints for the station
        """
        return self._constraints


    @property
    def bands(self) -> list[str]:
        """Observing bands the station can observe.
        Returns a dict_keys object with all bands in a string format as introduced in the freqs_sefd
        attribute when the Station was created.
        """
        return list(self._freqs_sefds.keys())


    @property
    def sefds(self) -> dict[str, u.Quantity]:
        """Returns a dictionary with the system equivalent flux density (SEFDs) for each
        of the frequencies the station can observe (given as keys).
        """
        return self._freqs_sefds


    @property
    def mount(self) -> Mount:
        """Returns the mount of the station
        """
        return self._mount


    @property
    def max_datarate(self) -> Union[u.Quantity, dict[str, u.Quantity], None]:
        """Returns the maximum data rate that the station can observe.
        It can be either a astropy.unit.Quantity value (equivalent to Mb/s),
        or a dictionary in the case that the antenna can show different data rates
        when it is participating in different networks. In that case, the abbreviation of the network will
        be the key of the dictionary, with the quantity as value. For example, Darnhall can record at 4 Gbps
        when participating within e-MERLIN observations, but their data rate is limited to 512 Mbps when
        participating within the EVN.
        """
        return self._max_datarate


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


    def altaz(self, obs_times: Time, target: FixedTarget) -> coord.SkyCoord:
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


    def hour_angle(self, obs_times: Time, target: FixedTarget) -> coord.Angle:
        """Returns the local hour angle of the target at the given time.
        """
        return self.observer.target_hour_angle(time=obs_times, target=target)


    def is_observable(self, obs_times: Time, target: FixedTarget) -> list[bool]:
        """Returns when the target source is observable for this station at the given times.

        Inputs
        - obs_times : astropy.time.Time
            Time to compute the elevation of the source
            (either a single timestamp or an array of times).
        - target : astroplan.FixedTarget
             Target coordinates to observe. If None, the target would be assumed to be visible at all times.

        Output
        - visible : list[bool]
            List of booleans of same length as targets for whether or not each target is ever
            observable in the time range given the constraints.
        """
        return is_observable(self.constraints, self.observer, target, times=obs_times)


    def is_always_observable(self, obs_times: Time, target: FixedTarget) -> bool:
        """Returns whenever the target source is observable for this station at all given times.

        Inputs
        - obs_times : astropy.time.Time
            Time to compute the elevation of the source
            (either a single timestamp or an array of times).
        - target : astroplan.FixedTarget
             Target coordinates to observe. If None, the target would be assumed to be visible at all times.

        Output
        - visible : bool
            Whether or not each target is observable in the time range given the constraints.
        """
        return is_always_observable(self.constraints, self.observer, target, times=obs_times)[0]


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


    def sefd(self, band: str) -> u.Quantity:
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
                 freqs_sefds: dict[str, u.Quantity],
                 fullname: Optional[str] = None, all_networks: Optional[list[str]] = None, country: str = '',
                 diameter: str = '', real_time: bool = False, mount: Optional[Mount] = None,
                 max_datarate: u.Quantity | dict[str, u.Quantity]| None = None,
                 selected: bool = True) -> None:
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
        - fullname : str [OPTIONAL]
            Full name of the station. If not given, same as `name` is assumed.
            It can be used to expand the full name if an abbreviation is typically used for the name.
            For example, name: VLA, fullname: Karl G. Jansky Very Large Array.
        - all_networks : list[str] [OPTIONAL]
            All networks where the station can participate.
            By default it will always include 'network' too.
        - country : str [OPTIONAL]
            Country where the station is located.
        - diameter : str [OPTIONAL]
            Diameter of the station (free format string). We recommend a syntax of e.g. '30 m' for normal
            single-dish antennas, and in case of interferometers it can have a form like '25 x 20 m',
            meaning that the station is composed of 25 antennas of 20 m each.
        - real_time : bool [OPTIONAL]
            If the station can participate in real-time observations (e.g. e-EVN), False by default.
        - mount : Mount [OPTIONAL]
            The mount of the station, including the type of mount and the slewing limits, speed, and acceleration.
            If not provided, it will assume an ALTAZ mount with no pointing limits and very high slewing speed.
        - max_datarate : u.Quantity or dict[str, u.Quantity]  [OPTIONAL]
            Specifies the maximum data rate that the station can record. If not provided, it will use the
            data rate assumed in the observation. It can be either a astropy.unit.Quantity value
            (equivalent to Mb/s), or a dictionary in the case that the antenna can show different data rates
            when it is participating in different networks. In that case, the abbreviation of the network will
            be the key of the dictionary, with the quantity as value. For example, Darnhall can record at 4 Gbps
            when participating within e-MERLIN observations, but their data rate is limited to 512 Mbps when
            participating within the EVN.
        - selected : bool [OPTIONAL]
            If the station is selected to participate in a given observation or not. True by default.
        """
        assert isinstance(selected, bool), "'selected' must be a bool"
        self._selected = selected
        super().__init__(name, codename, network, location,
                         freqs_sefds, fullname,
                         all_networks, country, diameter, real_time, mount, max_datarate)


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



class Network(object):
    """Defines a network (collection) of stations (`Station` objects) that
    can participate in an observation together.
    """
    def __init__(self, name: str, full_name: Optional[str] = None, stations: Optional[Iterable[Station]] = None,
                 observing_bands: Optional[Sequence[str]] = None,
                 max_datarates: Optional[Union[Sequence[u.Quantity], u.Quantity]] = None) -> None:
        """Initializes a Network of antennas.

        Inputs
        - name : str
            Name associated to the network of stations.
        - full_name : str  [OPTIONAL]
            Full name (expanded) associated to the network.
            If not provided, it will be name
        - stations : list of Station-type elements [OPTIONAL]
            List with all stations belonging to the given network.
            If not provided, it will initialize a network with zero antennas.
        - observing_bands : list of str  [OPTIONAL]
            List with all possible observing bands, as specified in the freqsetups.py file.
            If not provided, it will consider all observing bands from the available antennas.
        - max_datarates :  list of u.Quantity  or u.Quantity
            Maximum data rates to be able to record at each observing band. If provided as list,
            it must have the same dimensions as `observing_bands`. If scalar, means that the
            network has the same maximum data rate at all bands.
            If not provided, it will assume no upper limit on the data rate.
        """
        assert isinstance(name, str), "'name' must be a str."
        assert isinstance(stations, abc.Iterable) or stations is None, "'stations' must be a list or be empty."
        if observing_bands is None:
            assert not np.iterable(max_datarates), \
                   "`max_datarates` cannot be iterable if `observing_bands` is None."

        self._name = name
        self._fullname = name if full_name is None else full_name
        self._stations: dict[str, Station] = {}

        self._bands: dict[str, u.Quantity] = {}
        if (observing_bands is not None):
            # if np.iterable(max_datarates):
            if isinstance(max_datarates, Sequence):
                assert len(observing_bands) == len(max_datarates), \
                       "'observing_bands' and 'max_datarates' must have the same dimensions."
                for i,band in enumerate(observing_bands):
                    self._bands[band] = max_datarates[i]
            else:
                for i,band in enumerate(observing_bands):
                    self._bands[band] = max_datarates

        if stations is not None:
            for a_station in stations:
                assert isinstance(a_station, Station), \
                    f"There is an element in 'stations' that is not a Station object ({a_station})."
                if a_station.codename not in self._stations.keys():
                    self._stations[a_station.codename] = a_station
                else:
                    print(f"WARNING: {a_station.codename} is duplicated in the 'stations' list.")

                if observing_bands is None:
                    for aband in a_station.bands:
                        self._bands[aband] = max_datarates


    @property
    def name(self) -> str:
        """Name of the network of stations.
        """
        return self._name


    @name.setter
    def name(self, new_name: str):
        self._name = new_name


    @property
    def full_name(self) -> str:
        """Returns the full (expanded) name of the network.
        """
        return self._fullname

    @full_name.setter
    def full_name(self, new_full_name: str):
        self._fullname = new_full_name


    @property
    def stations(self) -> list[Station]:
        """Returns a list containing all stations in the network.
        """
        return list(self._stations.values())

    @stations.setter
    def stations(self, new_stations: list[Station]):
        self._stations = {s.codename: s for s in new_stations}

    @property
    def number_of_stations(self) -> int:
        """Returns the total number of stations in the network.
        """
        return len(self.stations)

    @property
    def names(self) -> list[str]:
        """Returns the names from all the stations in the network
        """
        return [s.name for s in self._stations.values()]

    @property
    def codenames(self) -> list[str]:
        """Returns a dict_keys with the `codenames` from all the stations in the network.
        """
        return list(self._stations.keys())


    @property
    def observing_bands(self) -> list[str]:
        """Returns a set with all `bands` that the stations in the network can observe,
        or at least a subset of stations.
        """
        return list(self._bands.keys())


    def max_datarate(self, observing_band: str) -> u.Quantity:
        return self._bands[observing_band]


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


    def __str__(self):
        return f"<{self.name}: <{self.number_of_stations}><{', '.join(self.codenames)}>>"


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


    @classmethod
    def get_stations_from_configfile(cls, filename: Optional[str] = None, codenames: Optional[list[str]] = None,
                                     name: str = 'network') -> Self:
        """Creates a Network object (i.e. a network of stations) by reading the station
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
            - mount = ALTAZ/EQUAT - type of the station mount.
            - ax1lim = (x0, x1) - Limits in degrees for the axis 1 of the mount.
            - ax2lim = (x0, x1) - Limits in degrees for the axis 2 of the mount.
            - ax1rate - Slewing speed in deg/s for the axis 1 of the mount.
            - ax2rate - Slewing speed in deg/s for the axis 2 of the mount.
            - ax1acc - Slewing acceleration in deg/s/s for the axis 1 of the mount.
            - ax2acc - Slewing acceleration in deg/s/s for the axis 2 of the mount.
            - real_time = yes/no - if the station can participate in real-time observations (e.g. e-EVN).
            - SEFD_**  - SEFD (in Jy units) of the station at the **cm band. If a given band is not present,
                        it is assumed that the station cannot observe it.
                        For example SEFD_21 = 500 means that the SEFD at 21cm is 500 Jy.
            - maxdatarate_** - maximum data rate of ** (in Mbps) that the station can achieve for the
              specify network name (e.g. EVN).
            - Any other attribute is accepted, but ignored in this code. That would easily allow future
              extensions of the code.
        - codenames : list
            If you only want to select a subset of all stations available in the input file,
            here you can pass a list with the codenames of the stations that should be imported.
        - name : str
            Name to assign to the network of stations that will be created.

        Returns
        - network : Network
            Returns a Network object containing the selected stations.
        """
        config = configparser.ConfigParser()
        if filename is None:
#             from importlib.resources import files
# data_text = files('mypkg.data').joinpath('data1.txt').read_text()
            with resources.as_file(resources.files("vlbiplanobs.data").joinpath("stations_catalog.inp")) \
                                                                     as stations_catalog_path:
                config.read(stations_catalog_path)
        else:
            # With this approach it raises a FileNotFound exception.
            # Otherwise config will run smoothly and provide an empty list.
            config.read(open(filename, 'r'))

        networks = cls(name, None, [])
        for stationname in config.sections():
            if (codenames is None) or (config[stationname]['code'] in codenames):
                temp = [float(i.strip()) for i in config[stationname]['position'].split(',')]
                a_loc = coord.EarthLocation(u.Quantity(temp[0], u.m), u.Quantity(temp[1], u.m),
                                            u.Quantity(temp[2], u.m))
                # Getting the SEFD values for the bands
                does_real_time = True if config[stationname]['real_time']=='yes' else False
                sefds = {}
                max_dt: u.Quantity | dict[str, u.Quantity] | None = None
                for akey in config[stationname].keys():
                    if 'SEFD_' in akey.upper():
                        sefds[f"{akey.upper().replace('SEFD_', '').strip()}"] = \
                                            float(config[stationname][akey])*u.Jy
                    if 'maxdatarate' == akey.lower():
                        max_dt = int(config[stationname][akey])*u.Mb/u.s
                    elif 'maxdatarate_' in akey.lower():
                        if not isinstance(max_dt, dict):
                            max_dt = {}

                        val = int(akey.removeprefix('maxdatarate_').strip())*u.Mb/u.s
                        net = [n.strip() for n in config[stationname][akey].split(',')]
                        for n in net:
                            max_dt[n] = val

                if config[stationname]['code'] in networks.codenames:
                    raise ValueError(f"The antenna with code {config[stationname]['code']} is "
                                     "duplicated in the input file.")

                if all([key in config[stationname] for key in ('mount', 'ax1rate', 'ax2rate', 'ax1lim',
                                                               'ax2lim')]):
                    configs = {}
                    try:
                        for alim in ('ax1lim', 'ax2lim'):
                            # Because they are the cable limits...
                            lims = [float(i) for i in config[stationname][alim].split(',')]
                            # easy check, if the range is >360, they can observe all azimuths
                            if (lims[1] - lims[0]) > 360:
                                configs[alim] = tuple((-1*u.deg, 361*u.deg))
                            else:
                                configs[alim] = tuple((lims[0]*u.deg, lims[1]*u.deg))

                        for arate in ('ax1rate', 'ax2rate'):
                            configs[arate] = u.Quantity(float(config[stationname][arate].strip()), u.deg/u.s)

                        if all([key in config[stationname] for key in ('ax1acc', 'ax2acc')]):
                            for aacc in ('ax1acc', 'ax2acc'):
                                configs[aacc] = u.Quantity(float(config[stationname][aacc].strip()), u.deg/u.s/u.s)

                            amount = Mount(MountType[config[stationname]['mount']],
                                           Axis(configs['ax1lim'], configs['ax1rate'], configs['ax1acc']),
                                           Axis(configs['ax2lim'], configs['ax2rate'],
                                                config[stationname]['ax2acc']))
                        else:
                            amount = Mount(MountType[config[stationname]['mount']],
                                           Axis(configs['ax1lim'], config[stationname]['ax1rate']),
                                           Axis(configs['ax2lim'], config[stationname]['ax2rate']))
                    except ValueError:
                        raise ValueError(f"when loading the data from antenna {stationname}.")

                else:
                    amount = None

                if (config[stationname]['possible_networks'] is not None) and \
                           (len(config[stationname]['possible_networks']) > 0):
                    the_networks: list[str] | None = [s.strip() \
                            for s in config[stationname]['possible_networks'].split(',')]
                else:
                    the_networks = None

                new_station = SelectedStation(stationname, config[stationname]['code'],
                        config[stationname]['network'], a_loc, sefds,
                        config[stationname]['station'], the_networks, config[stationname]['country'],
                        config[stationname]['diameter'], does_real_time, amount, max_dt)
                networks.add(new_station)

        if codenames is not None:
            for a_codename in codenames:
                if a_codename not in networks:
                    rprint(f"\n[yellow]WARNING: The antenna {a_codename} was not found in the catalogs.[/yellow]\n")
        return networks


    def stations_with_band(self, band: str, output_network_name: Optional[str] = None) -> Network:
        """Given the current network, it creates a sub-network including only the stations
        that can observe at the given band.

        Inputs
        - band : str
            The observing band that will be selected in the stations.
        - output_network_name : str [OPTIONAL]
            The name assigned to the new network.
            By default, it will be the same as the original network followed by @{band}.

        Returns
        - subnetwork : Network
            A subset of the original Network containing only the stations that can observe
            at the given band.
        """
        if output_network_name is None:
            output_network_name = f"{self.name}@{band}"

        subnetwork = Network(output_network_name, None, [])
        for station in self.stations:
            if band in station.bands:
                subnetwork.add(station)

        return subnetwork


    def select_stations(self, name: Optional[str] = None, codenames: Optional[list[str]] = None,
                        networknames: Optional[Union[str, list[str]]] = None,
                        full_name: Optional[str] = None) -> Network:
        """Returns a new Network object which will only contain the stations
        defined by the given list of codenames. It will thus be a subset of the current
        network.

        Input
        - name : str [OPTIONAL]
            Name assigned to the new network. If not provided, it will be the original
            name with the 'sub' preffix.
        - codenames : list[str]    [OPTIONAL]
            List with the codenames of the stations that should be present in the new
            network.
        - networknames : str | list[str]    [OPTIONAL]
            Network name or list of network names that should be present in the new network.
        - full_name : str [OPTIONAL]
            Full (expanded) name assigned to the new network. If not provided, it will use 'name'.
        Returns
        - subnetwork : Network
            A new Network object containing only the defined stations.
        Exceptions
        - It may raise KeyError if one of the given codenames are not present
          among the current stations.
        """
        subnetwork = Network(name if name is not None else f"sub{self.name}", full_name, [])

        if networknames is not None:
            pass

        elif codenames is not None:
            for codename in codenames:
                subnetwork.add(self[codename])


        return subnetwork


    @staticmethod
    def get_network_names_from_configfile(filename: Optional[str] = None) -> dict[str, Network]:
        """Reads a config file containing the different VLBI networks defined as a config parser file.
        Returns a dictionary with the nickname of the VLBI network as keys, and the information as
        keys in a second-order dict.

        Inputs
        - filename : str  OPTIONAL
            Path to the text file containing the information from all stations.
            If not provided, it reads the default station catalog file located in
                            data/network_catalog.inp
            Any other file should contain the same format (standard Python input config files),
            with the following fields per network (whose name would be provided as the name of
            section).
            - name - full name of the network.
            - default_antennas - comma-separated list with the codename of the antennas involved in the network.
            - max_datarate - integer number with the maximum datarate allowed in the network.
            - observing_bands - comma-separated list with the bands that can be observed (following
              the definition done in vlbiplanobs/freqsetups.py; e.g. 21, 18, etc, in cm).

        Returns
        - networks : dict[str, Network]
            Returns a dictionary containing the differnet networks.
        """
        config = configparser.ConfigParser()
        if filename is None:
            with resources.as_file(resources.files("vlbiplanobs.data").joinpath("network_catalog.inp")) \
                                                                     as networks_catalog_path:
                config.read(networks_catalog_path)
        else:
            # With this approach it raises a FileNotFound exception.
            # Otherwise config will run smoothly and provide an empty list.
            config.read(open(filename, 'r'))

        all_ants = Network.get_stations_from_configfile()
        networks: dict[str, Network] = dict()
        for networkname in config.sections():
            temp: str = config[networkname]['max_datarate']
            if ',' in temp:
                max_dt: list[u.Quantity] | u.Quantity = [int(dt.strip())*u.Mb/u.s for dt in temp.split(',')]
            else:
                max_dt = int(temp)*u.Mb/u.s

            obs_bands = [b.strip() for b in config[networkname]['observing_bands'].split(',')]
            default_ant = [a.strip() for a in config[networkname]['default_antennas'].split(',')]
            assert all([ant in all_ants for ant in default_ant]), "The default antenna(s) " \
                    f"({' ,'.join([ant for ant in default_ant if ant not in all_ants])}) from '{networkname}' " \
                    "is not present in stations_catalog!"

            assert all([band in freqsetups.bands.keys() for band in obs_bands]), \
                   f"Observing band ({', '.join([b for b in obs_bands if b not in freqsetups.bands.keys()])}) " \
                   "not present in freqsetups.py!"

            if isinstance(max_dt, Sequence):
                assert all([int(dr.value) in freqsetups.data_rates.keys() for dr in max_dt]), \
                       f"Data rate ({', '.join([d for d in max_dt if int(d.value) not in freqsetups.data_rates.keys()])}) not present in freqsetups.py!"
            else:
                assert int(max_dt.value) in freqsetups.data_rates.keys(), f"Data rate ({max_dt}) not " \
                       "present in freqsetups.py!"

            antennas = [ant for ant in all_ants if networkname in ant.all_networks]
            # antennas = [all_ants[ant] for ant in default_ant]
            for ant in antennas:
                ant.selected = ant.codename in default_ant

            assert len(antennas) > 0, f"No antennas found for the network {networkname}."

            networks[networkname] =  Network(name=networkname, full_name=config[networkname]['name'],
                                             stations=antennas,
                                             observing_bands=obs_bands,
                                             max_datarates=max_dt)

        return networks

