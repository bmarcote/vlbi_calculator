

class Station(object):

    def __init__(self, name, codename, network, location, freqs_sefds, min_elevation=20*u.deg,
                 fullname=None, all_networks=None, country='', diameter='', real_time=False):
        """Initializes a station. The given name must be the name of the station that
        observes, with the typical 2-letter format used in the EVN (with exceptions).

        Inputs
        - name : str
            Name of the observer (the station that is going to observe).
        - codename : str
            A code for the name of the station. It can be the same as name.
        - network : str
            Name of the network to which the station belongs.
        - location : EarthLocation
            Position of the observer on Earth.
        - freqs_sefds : dict
            Dictionary with all frequencies the station can observe, and as values
            the SEFD at each frequency.
        - min_elevation : Quantity
            Minimum elevation that the station can observe a source. If no units
            provided, degrees are assumed. By default it 20 degrees.
        - fullname : str [OPTIONAL]
            Full name of the station. If not given, `name` is assumed.
        - all_networks : str [OPTIONAL]
            Networks where the station can participate (free style).
        - country : str [OPTIONAL]
            Country where the station is placed.
        - diameter : str [OPTIONAL]
            Diameter of the station (free format).
        - real_time : bool [OPTIONAL, False by default]
            If the station can participate in real-time observations (e.g. e-EVN).
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
    def name(self):
        """Name of the station.
        """
        return self.observer.name

    @property
    def codename(self):
        """Codename of the station (typically a two-letter accronym).
        """
        return self._codename

    @property
    def fullname(self):
        return self._fullname

    @property
    def network(self):
        """Name of the network to which the station belongs.
        """
        return self._network

    @property
    def all_networks(self):
        """Name of all networks to which the station belongs.
        """
        return self._all_networks

    @property
    def country(self):
        return self._country

    @property
    def diameter(self):
        return self._diameter

    @property
    def real_time(self):
        return self._real_time

    @property
    def location(self):
        """Location of the station in EarthLocation type.
        """
        return self.observer.location

    @property
    def bands(self):
        """Observing bands the station can observe.
        """
        return self._freqs_sefds.keys()

    @property
    def sefds(self):
        """Returns a dictionary with the SEFDs for each of the frequencies
        the station can observe (given as keys).
        """
        return self._freqs_sefds

    @property
    def min_elevation(self):
        """Minimum elevation the station can observe a source.
        """
        return self._min_elev


    def elevation(self, obs_times, target):
        """Returns the elevation of the source as seen by the Station during obs_times.

        Inputs
        - obs_times : astropy.time.Time
            Time to compute the elevation of the source (either single time or a list of times).
        - target : astroplan.FixedTarget
             Target to observe.

        Output
        - elevations : ndarray
            Elevation of the source at the given obs_times
        """
        # source_altaz = source_coord.transform_to(coord.AltAz(obstime=obs_times,
        #                                                 location=self.location))
        # return source_altaz.alt
        return self.observer.altaz(obs_times, target).alt


    def altaz(self, obs_times, target):
        """Returns the altaz coordinates of the target for the given observing times.
        """
        return self.observer.altaz(obs_times, target)


    def is_visible(self, obs_times, target):
        """Return if the source is visible for this station at the given times.
        """
        elevations = self.elevation(obs_times, target)
        return np.where(elevations >= self.min_elevation)

    def has_band(self, band):
        """Returns if the Station can observed the given band `the_band`.
        """
        return band in self.bands

    def sefd(self, band):
        """Returns the SEFD of the Station at the given band.
        """
        return self._freqs_sefds[band]

    def __str__(self):
        return f"<{self.codename}>"

    def __repr__(self):
        return f"<stations.Station: {self.codename}>"



class SelectedStation(Station):
    def __init__(self, name, codename, network, location, freqs_sefds, min_elevation=20*u.deg,
                 fullname=None, all_networks=None, country='', diameter='', real_time=False, selected=True):
        self._selected = selected
        super().__init__(name, codename, network, location, freqs_sefds,
                         min_elevation, fullname, all_networks, country, diameter, real_time)

    @property
    def selected(self):
        return self._selected

    @selected.setter
    def selected(self, isselected):
        assert isinstance(isselected, bool)
        self._selected = isselected


class Stations(object):
    """Class that collects groups of stations (`Station` objects) that form
    an observatory/array/network.
    """
    def __init__(self, name, stations):
        """Inputs
            - name : str
                Name associated to the observatory/array/network.
            - stations : list of Stations
                List with all stations belonging to the network.
        """
        self._name = name
        self._stations = {}
        for a_station in stations:
            assert a_station.codename not in self._stations.keys()
            self._stations[a_station.codename] = a_station

        self._keys = tuple(self._stations.keys())

    @property
    def name(self):
        """Name of the observatory/array/network.
        """
        return self._name

    @property
    def stations(self):
        """Returns the list of all stations in the observatory/array/network.
        """
        return list(self._stations.values())

    @property
    def number_of_stations(self):
        return len(self.stations)

    def add(self, a_station):
        """Adds a new station to the observatory/array/network.
        """
        if a_station.codename in self._stations.keys():
            print(f"WARNING: {a_station.codename} already in {self.name}.")

        self._stations[a_station.codename] = a_station
        self._keys = tuple(self._stations.keys())

    def keys(self):
        """Returns a tuple of `codenames` from the stations.
        """
        return self._keys

    def __str__(self):
        return f"<{self.name}: <{', '.join(self._stations.keys())}>>"

    def __len__(self):
        return self._stations.__len__()

    def __getitem__(self, key):
        if isinstance(key, int):
            return self._stations[self.keys()[key]]
        else:
            return self._stations[key]

    def __setitem__(self, key, value):
        self._stations[key] = value
        self._keys = tuple(self._stations.keys())

    def __delitem__(self, key):
        if isinstance(key, int):
            self._stations.__delitem__(self.keys[key])
        else:
            self._stations.__delitem__(key)

        self._keys = tuple(self._stations.keys())

    def __iter__(self):
        return iter(self._stations.values())

    def __contains__(self, item):
        return self._stations.__contains__(item)





