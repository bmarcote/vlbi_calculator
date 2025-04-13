# -*- coding: utf-8 -*-
# Licensed under GPLv3+ - see LICENSE
from __future__ import annotations
from collections import abc
from typing import Optional, Union, Iterable, Sequence, Self, Generator
import configparser
from importlib import resources
import numpy as np
from rich import print as rprint
from astropy import units as u
from astropy import coordinates as coord
from astropy.time import Time
from astroplan import Observer, FixedTarget, is_observable, is_event_observable, is_always_observable
from dataclasses import dataclass
from enum import Enum, auto
from . import constraints
from . import freqsetups

"""Module that defines the `Station` and `Stations` objects, which represent a station (antenna)
or a network composed of antennas.
"""

__all__: list[str] = ["Station", "Stations", "Mount", "MountType"]


class MountType(Enum):
    ALTAZ = auto()
    EQUAT = auto()


@dataclass
class Axis:
    """Defines one of the mount axis for a given antenna.
    The axis is defined by the angle range it can move (limits), and the slewing speed and
    acceleration.
    """
    limits: tuple[u.Quantity, u.Quantity]
    speed: u.Quantity
    acceleration: u.Quantity = u.Quantity(0.0, u.deg / u.s / u.s)


@dataclass
class Mount:
    """Defines a telescope mount, including the mount type (ALTAZ, EQUAT), the limits of the
    telescope to point in the sky, and the slewing velocity (rate) and acceleration (acc).
    """
    mount_type: MountType
    ax1: Axis
    ax2: Axis


class Station(object):
    """Defines an astronomical station (antenna).
    A station is defined by some names and codenames, coordinates, and its sensitivity for the
    radio bands (wavelengths) it can observe.
    Apart of the metadata related to the station, it allows to compute the altitude/azimuth,
    elevation, or simply when a source is visible from the station for a given time range.
    """

    def __init__(self, name: str, codename: str, networks: tuple[str], location: coord.EarthLocation,
                 freqs_sefds: dict[str, u.Quantity], fullname: Optional[str] = None,
                 country: str = "", diameter: str = "",
                 real_time: bool = False, mount: Optional[Mount] = None,
                 max_datarate: Optional[u.Quantity | dict[str, u.Quantity]] = None,
                 datarate: Optional[u.Quantity] = None) -> None:
        """Initializes a station.

        Inputs
        - name : str
            Name of the station.
            If it contains undercores (_), they will be converted to blank spaces.
        - codename : str
            A short code (accronym) for the name of the station. It is meant to follow the standard
            approach from the EVN: an (often) two-letter code unique for each station.
        - networks : tuple[str]
            Name of the network(s) to which the station belongs (e.g. EVN).
        - location : astropy.coordinates.EarthLocation
            Position of the station on Earth in (x,y,z) gecentric coordinates.
        - freqs_sefds : dict
            Dictionary with all frequencies the station can observe as keys of the dictionary, and
            the values representing the system equivalent flux density (SEFD; in Jansky units)
            at each frequency.
            Although the key format is in principle free, we recommend to use the syntax 'XXcm' (str type).
            This will be then consistent with the default station catalog and will avoid issues for some
            functions.
        - fullname : str [OPTIONAL]
            Full name of the station. If not given, same as `name` is assumed.
            It can be used to expand the full name if an abbreviation is typically used for the name.
            For example, name: VLA, fullname: Karl G. Jansky Very Large Array.
        - country : str [OPTIONAL]
            Country where the station is located.
        - diameter : str [OPTIONAL]
            Diameter of the station (free format string). We recommend a syntax of e.g. '30 m' for normal
            single-dish antennas, and in case of interferometers it can have a form like '25 x 20 m',
            meaning that the station is composed of 25 antennas of 20 m each.
        - real_time : bool [OPTIONAL]
            If the station can participate in real-time observations (e.g. e-EVN), False by default.
        - mount : Mount [OPTIONAL]
            The mount of the station, including the type of mount and the slewing limits, speed,
            and acceleration.
            If not provided, it will assume an ALTAZ mount with no pointing limits and very high
            slewing speed.
        - max_datarate : u.Quantity or dict[str, u.Quantity]  [OPTIONAL]
            Specifies the maximum data rate that the station can record. If not provided, it will use the
            data rate assumed in the observation. It can be either a astropy.unit.Quantity value
            (equivalent to Mb/s), or a dictionary in the case that the antenna can show different data rates
            when it is participating in different networks. In that case, the abbreviation of the
            network will be the key of the dictionary, with the quantity as value.
            For example, Darnhall can record at 4 Gbps when participating within e-MERLIN observations,
            but their data rate is limited to 512 Mbps when participating within the EVN.
        - datarate : u.Quantity [OPTIONAL]
            Specifies the data rate that the station will record. If not provided, it will use the
            data rate assumed in the observation.
        """
        for a_var, a_var_name in zip((name, codename, country, diameter),
                                     ("name", "codename", "country", "diameter")):
            assert isinstance(a_var, str), f"'{a_var_name}' must be a str."

        assert isinstance(fullname, str) or fullname is None
        assert isinstance(real_time, bool), "'real_time' must be a bool."
        self.observer: Observer = Observer(name=name.replace("_", " "), location=location)
        self._codename: str = codename
        self._networks: tuple[str] = networks if isinstance(networks, tuple) else (networks,)
        self._freqs_sefds: dict[str, float] = {f if "cm" in f else f"{f}cm": v
                                               for f, v in freqs_sefds.items()}
        self._fullname: str = name if fullname is None else fullname
        self._country: str = country
        self._diameter: str = diameter
        self._real_time: bool = real_time
        if mount is None:
            self._mount: Mount = Mount(MountType.ALTAZ, Axis((-10*u.deg, 370*u.deg),
                                                             300*u.deg/u.s, 0.0*u.deg/u.s/u.s),
                                       Axis((10*u.deg, 100*u.deg), 300*u.deg/u.s, 0.0*u.deg/u.s/u.s))
        else:
            assert isinstance(mount, Mount)
            self._mount = mount

        self._max_datarate: u.Quantity | dict | None = max_datarate
        self._datarate: u.Quantity | None = datarate
        if self.mount.mount_type == MountType.ALTAZ:
            self._constraints = [constraints.AzimuthConstraint(min=self.mount.ax1.limits[0],
                                                               max=self.mount.ax1.limits[1]),
                                 constraints.ElevationConstraint(min=self.mount.ax2.limits[0],
                                                                 max=self.mount.ax2.limits[1])]
        else:
            self._constraints = [constraints.HourAngleConstraint(min=self.mount.ax1.limits[0],
                                                                 max=self.mount.ax1.limits[1]),
                                 constraints.DeclinationConstraint(min=self.mount.ax2.limits[0],
                                                                   max=self.mount.ax2.limits[1]),
                                 constraints.ElevationConstraint(min=5*u.deg)]

    @property
    def name(self) -> str:
        """Name of the station."""
        assert self.observer.name is not None
        return self.observer.name

    @property
    def codename(self) -> str:
        """Codename of the station (typically a two-letter accronym)."""
        return self._codename

    @property
    def fullname(self) -> str:
        """Full name of the station. If not specified, it can be the same as 'name'."""
        return self._fullname

    @property
    def networks(self) -> tuple[str]:
        """Name of the network(s) to which the station belongs."""
        return self._networks

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
        """If the station can participate in real-time observations (e.g. e-EVN)."""
        return self._real_time

    @property
    def location(self) -> coord.EarthLocation:
        """Location of the station as an astropy.coordinates.EarthLocation object."""
        return self.observer.location

    @property
    def constraints(self) -> list[constraints.Constraint]:
        """Returns the observing constraints for the station"""
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
        """Returns the mount of the station"""
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

    @property
    def datarate(self) -> u.Quantity | None:
        """Returns the data rate that the station can observe.
        It can be either a astropy.unit.Quantity value (equivalent to Mb/s) or None (if not specified).
        """
        return self._datarate

    @datarate.setter
    def datarate(self, new_datarate: u.Quantity):
        self._datarate = new_datarate

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
        """Returns the local hour angle of the target at the given time."""
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
        - visible : list[list[bool]]
            List of booleans of same length as targets for whether or not each target is ever
            observable in the time range given the constraints.
        """
        return is_event_observable(self.constraints, self.observer, target, times=obs_times)[0, :]

    def is_ever_observable(self, obs_times: Time, target: FixedTarget) -> list[bool]:
        """Returns if the target source is observable for this station at least for part of the times.

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

    def slewing_time(self, coordinates1: coord.AltAz, coordinates2: coord.AltAz, time: Time) -> u.Quantity:
        # TODO: this is not the best implementation (asking for a time). Check later.
        """Returns the expected slewing time to move from coordinates1 to coordinates2 according to
        the known velocity and acceleration from the mount.

        Args
        - coordinates1 : astropy.coordinates.AltAz
            AltAz object with the coordinates of the first point. Must be in altaz coordinates.
        - coordinates2 : astropy.coordinates.AltAz
            AltAz object with the coordinates of the second point. Must be in altaz coordinates.

        Returns
        - slewing_time : astropy.time.TimeDelta
            TimeDelta object with the slewing time.
        """
        assert isinstance(coordinates1, coord.AltAz) and isinstance(coordinates2, coord.AltAz)
        # If the antenna mount is AltAz, then all the same; but if it's equatorial, then it needs to be RA/DEC.
        separation = coordinates1.spherical_offset_to(coordinates2)
        if self.mount.mount_type == MountType.ALTAZ:
            altaz1 = self.altaz(time, coordinates1)
            altaz2 = self.altaz(time, coordinates2)
            separation = (altaz1.az - altaz2.az, altaz1.alt - altaz2.alt)
        else:
            separation = coordinates1.spherical_offset_to(coordinates2)

        t = [0, 0]
        for i, axis in enumerate([self.mount.ax1, self.mount.ax2]):
            if axis.acceleration == 0.0:
                t[i] = separation[i]/axis.speed
            elif separation[i] <= 0.5*axis.speed**2/axis.acceleration:
                t[i] = np.sqrt(2*separation[i]/axis.acceleration)
            else:
                t[i] = separation[i]/axis.speed + 0.5*axis.speed/axis.acceleration

        return max(t)

    def __str__(self):
        return f"<{self.codename}>"

    def __repr__(self):
        return f"<Station: {self.codename}>"


class Stations(object):
    """Defines a network of antennas (collection) that can participate in an observation together.
    """

    def __init__(self, filename: Optional[str] = None, stations: Optional[Iterable[Station]] = None,
                 observing_bands: Optional[Sequence[str]] = None,
                 max_datarates: Optional[Union[Sequence[u.Quantity], u.Quantity]] = None,
                 name: Optional[str] = None) -> None:
        """Initializes a Stations of antennas.

        Inputs
        - filename : str  [OPTIONAL]
            Name of the file that contains the information of the antennas
            (by default, the 'stations_catalog.inp').
        - stations : list of Station-type elements  [OPTIONAL]
            List with all stations belonging to the given network. If not provided, it will include
            all antennas present in the catalogs.
        - observing_bands : list of str  [OPTIONAL]
            List with all possible observing bands, as specified in the freqsetups.py file.
            If not provided, it will consider all observing bands from the available antennas.
        - name : str  [OPTIONAL]
            Full name (expanded) associated to the network.
            If not provided, it will be name
        - max_datarates :  list of u.Quantity  or u.Quantity or None [OPTIONAL]
            Maximum data rates to be able to record at each observing band. If provided as list,
            it must have the same dimensions as `observing_bands`. If scalar, means that the
            network has the same maximum data rate at all bands.
            If not provided, it will assume no upper limit on the data rate.

        Raises
            ValueError: If 'max_datarates' is provided, 'observing_bands' must also be provided.
        """
        assert isinstance(stations, abc.Iterable) or stations is None, "'stations' must be a list or be empty."
        self._stations: dict[str, Station] = {}
        self._bands: dict[str, u.Quantity] = {}
        self._name: Optional[str] = name
        if (observing_bands is None) and (max_datarates is not None):
            raise ValueError("If 'max_datarates' is provided, 'observing_bands' must also be provided.")

        if observing_bands is not None:
            if isinstance(max_datarates, Sequence):
                assert len(observing_bands) == len(max_datarates), \
                       "'observing_bands' and 'max_datarates' must have the same dimensions."
                for i, band in enumerate(observing_bands):
                    self._bands[band] = max_datarates[i]
            elif max_datarates is not None:
                for band in observing_bands:
                    self._bands[band] = max_datarates
            else:
                for band in observing_bands:
                    self._bands[band] = None

        if stations is not None:
            for a_station in stations:
                assert isinstance(a_station, Station), \
                       f"There is an element in 'stations' that is not a Station object ({a_station})."
                if a_station.codename not in self._stations.keys():
                    self._stations[a_station.codename] = a_station
                else:
                    rprint(f"[yellow bold]WARNING: {a_station.codename} is duplicated in the provided "
                           "stations list while initializing a network.[/yellow bold]")

                if observing_bands is None:
                    for aband in a_station.bands:
                        self._bands[aband] = None
        else:
            # Read all stations from catalog,
            for station in self._get_stations_from_configfile(filename):
                self._stations[station.codename] = station
                for aband in station.bands:
                    self._bands[aband] = None
            # I don't think I need to know the networks at this point
            # _all_networks: dict[str, Stations] = self.get_networks_from_configfile(networks_filename)
            # TODO: implement this part!!!!!!!!!!!!
            # raise NotImplementedError

    @property
    def name(self) -> Optional[str]:
        return self._name

    @property
    def stations(self) -> list[Station]:
        """Returns a list containing all stations in the network."""
        return list(self._stations.values())

    @stations.setter
    def stations(self, new_stations: list[Station]):
        self._stations = {s.codename: s for s in new_stations}

    @property
    def number_of_stations(self) -> int:
        """Returns the total number of stations in the network."""
        return len(self.stations)

    @property
    def station_names(self) -> list[str]:
        """Returns the names from all the stations in the network"""
        return [s.name for s in self._stations.values()]

    @property
    def station_codenames(self) -> list[str]:
        """Returns a dict_keys with the `codenames` from all the stations in the network."""
        return list(self._stations.keys())

    @property
    def observing_bands(self) -> list[str]:
        """Returns a set with all `bands` that the stations in the network can observe,
        or at least a subset of stations.
        """
        return list(self._bands.keys())

    def max_datarate(self, observing_band: str) -> Optional[u.Quantity]:
        return self._bands[observing_band]

    def add_station(self, station: Union[Station, Iterable[Station]]):
        """Adds a new station to the network.
        If a station with the same codename is already present, it will do nothing.

        Inputs
        - station : Station
            Station to be added to the network.
        """
        assert isinstance(station, Station) or isinstance(station, Iterable), \
               "station must be a Station object or a Iterable of Station objects."
        if isinstance(station, Station):
            if station.codename in self.station_codenames:
                rprint(f"[yellow bold]WARNING: {station.codename} already in stations. "
                       "Ignoring addition.[/yellow bold]")
            else:
                self._stations[station.codename] = station
        else:
            for a_stat in station:
                if a_stat.codename in self.station_codenames:
                    rprint(f"[yellow bold]WARNING: {a_stat.codename} already in stations. "
                           "Ignoring addition.[/yellow bold]")
                else:
                    self._stations[a_stat.codename] = a_stat

    def remove_station(self, a_station: Station):
        """Removes a station from the network.
        If the station is not present, then it will do nothing.

        Inputs
        - a_station : Station
            Station to be added to the network.
        """
        if a_station.codename in self.station_codenames:
            self._stations.__delitem__(a_station.codename)

    def __str__(self):
        return f"<Stations ({self.number_of_stations}): <{', '.join(self.station_codenames)}>>"

    def __len__(self):
        return self._stations.__len__()

    def __getitem__(self, key):
        if isinstance(key, int):
            return self._stations[self.station_codenames[key]]
        elif key in self._stations.keys():
            return self._stations[key]
        elif key in self.station_names:
            return self._stations[self.station_codenames[self.station_names.index(key)]]
        else:
            raise KeyError(f"Key {key} is not in the network.")

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

    def __add__(self, stations: Stations):
        if not isinstance(stations, Stations):
            raise TypeError("The element to add must be a Stations object")

        obs_bands = [b for b in self.observing_bands if b in stations.observing_bands]
        max_dt = []
        for b in obs_bands:
            if self.max_datarate(b) is None:
                max_dt.append(stations.max_datarate(b))
            elif stations.max_datarate(b) is None:
                max_dt.append(self.max_datarate(b))
            else:
                max_dt.append(min(self.max_datarate(b), stations.max_datarate(b)))  # type: ignore

        return Stations(stations=set(self.stations + stations.stations),
                        observing_bands=obs_bands, max_datarates=max_dt)

    def __iadd__(self, stations: Stations):
        return self.__add__(stations)

    def stations_with_band(self, band: str) -> Generator[Station]:
        """Returns a the stations in the Stations that can observe at the given band.
        that can observe at the given band.

        Inputs
        - band : str
            The observing band.

        Returns
        - stations : Generator[Station]
            A list containing only the stations that can observe at the given band.
        """
        for station in self.stations:
            if band in station.bands:
                yield station

    def filter_networks(self, networks: Union[str, Sequence[str]], only_defaults: bool = False) -> Stations:
        """Returns a new network which only contains the antennas that are part
        of the given network

        Input
            - networks : str or Sequence[str]
                The name of the network or networks to which all stations need to belong to.
            - only_defaults : bool  [DEFAULT = False]
                If True, it will only return the antennas that are part of the default
                network, ignoring all others.

        Returns
            - Stations
                A new network object containing only the antennas that can observe in the network.
        """
        # TODO: change this to network -> networks so it may include multiple!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if isinstance(networks, str):
            networks = [n.strip() for n in networks.split(',')]

        if only_defaults:
            all_networks = self.get_networks_from_configfile()
            if len(networks) > 0:
                the_network = all_networks[networks[0]]
                for a_network in networks[1:]:
                    the_network += all_networks[a_network]
            else:
                raise ValueError("networks cannot be an empty sequence!")

            bands, datarates = None, None
            for aband, adatarate in zip(self.observing_bands, self._bands.values()):
                if aband in the_network.observing_bands:
                    if bands is None:
                        bands = [aband,]
                    else:
                        bands.append(aband)
                    if datarates is None:
                        datarates = [adatarate,]
                    else:
                        datarates.append(adatarate)

            return Stations(stations=(s for s in self.stations if any(n in s.networks for n in networks)
                                      and s.codename in the_network.station_codenames),
                            observing_bands=bands, max_datarates=datarates)

        return Stations(stations=(s for s in self.stations if any(n in s.networks for n in networks)),
                        observing_bands=self.observing_bands, max_datarates=self._bands.values())

    def filter_band(self, band: str) -> Stations:
        """Returns a new network which only contains the antennas that can observe
        at the given observing band.

        Input
            - band : str
                The observing band.

        Returns
            - Stations
                A new network object containing only the antennas that can observe at the given band.
        """
        return Stations(stations=(s for s in self.stations if band in s.bands), observing_bands=(band,),
                        max_datarates=self.max_datarate(band))

    def filter_antennas(self, codenames: Sequence[str]) -> Stations:
        """Returns a new network object which will only contain the stations
        defined by the given list of codenames. It will thus be a subset of the current
        network.

        Input
        - codenames : Sequence[str]
            List with the codenames of the stations that should be present in the new
            network.

        Returns
        - subnetwork : Stations
            A new Stations object containing only the defined stations.
        Exceptions
        - It may raise KeyError if one of the given codenames are not present
          among the current stations.
        """
        if not all([codename in self.station_codenames for codename in codenames]):
            unexpected_ant = set(codenames).difference(set(self.station_codenames))
            if len(unexpected_ant) == 1:
                rprint(f"[yellow bold]WARNING: The antenna with codename {unexpected_ant.pop()} is not present "
                       "in the current network.[/yellow bold]")
            else:
                rprint(f"[yellow bold]WARNING: The antennas with codenames {', '.join(list(unexpected_ant))} "
                       "are not present in the current network.[/yellow bold]")
        codename_indices = {codename: i for i, codename in enumerate(codenames)}
        return Stations(stations=sorted([s for s in self.stations if s.codename in codenames],
                                        key=lambda x: codename_indices[x.codename]),
                        observing_bands=self.observing_bands, max_datarates=self._bands.values())

    @staticmethod
    def get_networks_from_configfile(filename: Optional[str] = None) -> dict[str, Stations]:
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
        - networks : dict[str, Stations]
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
            config.read(open(filename, "r"))

        all_ants: Stations = Stations()
        networks: dict[str, Stations] = dict()
        for networkname in config.sections():
            temp: str = config[networkname]["max_datarate"]
            if "," in temp:
                max_dt: list[u.Quantity] | u.Quantity = [int(dt.strip())*u.Mb/u.s for dt in temp.split(",")]
            else:
                max_dt = int(temp)*u.Mb/u.s

            obs_bands = [b.strip() for b in config[networkname]["observing_bands"].split(",")]
            default_ant = [a.strip() for a in config[networkname]["default_antennas"].split(",")]
            assert all([ant in all_ants.station_codenames for ant in default_ant]), "The default antenna(s) " \
                   f"({', '.join([ant for ant in default_ant if ant not in all_ants])}) from '{networkname}' " \
                   "is not present in stations_catalog!"

            assert all([band in freqsetups.bands.keys() for band in obs_bands]), \
                   f"Observing band ({', '.join([b for b in obs_bands if b not in freqsetups.bands.keys()])}) " \
                   "not present in freqsetups.py!"

            if isinstance(max_dt, Sequence):
                temp = ', '.join([d for d in max_dt if int(d.value) not in freqsetups.data_rates.keys()])
                assert all([int(dr.value) in freqsetups.data_rates.keys() for dr in max_dt]), \
                       f"Data rate ({temp}) not present in freqsetups.py!"
            else:
                assert int(max_dt.value) in freqsetups.data_rates.keys(), \
                       f"Data rate ({max_dt}) not present in freqsetups.py!"

            antennas = [all_ants[[a.codename for a in all_ants].index(ant)] for ant in default_ant]
            assert len(antennas) > 0, f"No antennas found for the network {networkname}."

            networks[networkname] = Stations(stations=antennas, observing_bands=obs_bands, max_datarates=max_dt,
                                             name=config[networkname]["name"])

        return networks

    @staticmethod
    def _parse_station_from_configfile(stationname: str, station: dict) -> Station:
        """Parses the entry of a station from the config file and returns a Station object.

        Inputs
        - stationname : str
            The name of the station as defined in the header of the entry for
        - station : dict
            A dictionary containing the different information from a station, with the following
            fields per station (whose name would be provided as the name of section):
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

        Returns
            Station
        """
        temp = [float(i.strip()) for i in station["position"].split(",")]
        a_loc = coord.EarthLocation(u.Quantity(temp[0], u.m), u.Quantity(temp[1], u.m), u.Quantity(temp[2], u.m))
        # Getting the SEFD values for the bands
        does_real_time = True if station["real_time"] == "yes" else False
        sefds = {}
        max_dt: u.Quantity | dict[str, u.Quantity] | None = None
        for akey in station.keys():
            if "SEFD_" in akey.upper():
                sefds[f"{akey.upper().replace('SEFD_', '').strip()}"] = float(station[akey])*u.Jy
            if "maxdatarate" == akey.lower():
                max_dt = int(station[akey])*u.Mb/u.s
            elif "maxdatarate_" in akey.lower():
                if not isinstance(max_dt, dict):
                    max_dt = {}

                val = int(akey.removeprefix("maxdatarate_").strip())*u.Mb/u.s
                net = [n.strip() for n in station[akey].split(",")]
                for n in net:
                    max_dt[n] = val

        if all([key in station for key in ("mount", "ax1rate", "ax2rate", "ax1lim", "ax2lim")]):
            configs = {}
            try:
                for alim in ("ax1lim", "ax2lim"):
                    # Because they are the cable limits...
                    lims = [float(i) for i in station[alim].split(",")]
                    # easy check, if the range is >360, they can observe all azimuths
                    if (lims[1] - lims[0]) > 360:
                        configs[alim] = tuple((-1*u.deg, 361*u.deg))
                    else:
                        configs[alim] = tuple((lims[0]*u.deg, lims[1]*u.deg))

                if station["mount"] == 'EQUAT':
                    configs["ax1lim"] = tuple((ax1lim.value*u.hourangle for ax1lim in configs["ax1lim"]))

                for arate in ("ax1rate", "ax2rate"):
                    configs[arate] = u.Quantity(float(station[arate].strip()), u.deg/u.s)

                if all([key in station for key in ("ax1acc", "ax2acc")]):
                    for aacc in ("ax1acc", "ax2acc"):
                        configs[aacc] = u.Quantity(float(station[aacc].strip()), u.deg/u.s/u.s)

                    amount = Mount(MountType[station["mount"]], Axis(configs["ax1lim"],
                                                                     configs["ax1rate"],
                                                                     configs["ax1acc"]),
                                   Axis(configs["ax2lim"], configs["ax2rate"],
                                        station["ax2acc"]))
                else:
                    amount = Mount(MountType[station["mount"]], Axis(configs["ax1lim"],
                                   station["ax1rate"]), Axis(configs["ax2lim"],
                                   station["ax2rate"]))
            except ValueError:
                raise ValueError(f"when loading the data from antenna {station['station']}.")
        else:
            amount = None

        return Station(stationname, station["code"], tuple(station["networks"].split(",")),  # type: ignore
                       a_loc, sefds, station["station"], station["country"], station["diameter"],
                       does_real_time, amount, max_dt)

    @staticmethod
    def _get_stations_from_configfile(filename: Optional[str] = None,
                                      codenames: Optional[list[str]] = None) -> Generator[Station]:
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

        Returns
        - network : Generator[Station]
            Returns a Stations object containing the specified stations.
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
            config.read(open(filename, "r"))

        for stationname in config.sections():
            if (codenames is None) or (config[stationname]["code"] in codenames):
                yield Stations._parse_station_from_configfile(stationname, config[stationname])  # type: ignore
