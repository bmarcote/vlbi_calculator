import numpy as np
from astropy import units as u
from astropy import coordinates as coord
from astropy.time import Time, TimeDelta
from astroplan import FixedTarget

from vlbiplanobs.stations import Stations

__all__ = ['SourceNotVisible', 'Source', 'Observation']

"""Defines an observation, which basically consist of a given network of stations,
observing a target source for a given time range and at an observing band.
"""

class SourceNotVisible(Exception):
    """Exception produced when a given target source cannot be observed for any
    antenna in the network.
    """
    pass



class Source(FixedTarget):
    """Defines a target source located at some coordinates and with a given name.
    """
    def __init__(self, coordinates: str = None, name: str = None):
        """Initializes a Source object.

        Inputs
        - coordinates : str [OPTIONAL]
            Coordinates of the target source in a str format recognized by
            astropy.coordinates.SkyCoord (e.g. XXhXXmXXs XXdXXmXXs).
            J2000 coordinates are assumed.
            If not provided, name must be a source name recognized by astroquery.
        - name : str [OPTIONAL]
            Name associated to the source. By default is None.

        If both provided, the given coordinates will be used for the given source.

        It may raise:
        - NameResolveError: if the name is not recognized (and no coordinates are provided)
        - ValueError: if the coordinates have an unrecognized format.
        - AssertionError: if none, coordinates nor name, are provided.
        """
        assert (name is not None) or (coordinates is not None), \
               "At least one 'coordiantes' or 'name' must be provided"
        if (name is not None) and (coordinates is None):
            coordinates = coord.get_icrs_coordinates(name)

        super().__init__(coord.SkyCoord(coordinates), name)



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
    def __init__(self, target: Source = None, times: Time = None, band: str = None,
                 datarate=None, subbands: int = None, channels: int = None,
                 polarizations: int = None, inttime=None, ontarget: float = 1.0,
                 stations: Stations = None, bits: int = 2):
        """Initializes an observation.
        Note that you can initialize an empty observation at this stage and add the
        information for the different attributes later. However, you may raise exception
        if running methods that require some of the unset attributes.

        Inputs
        - target : vlbiplanobs.observation.Source
            Target source to be observed.
        - times : astropy.time.Time
            An array of times defining the duration of the observation. That is,
            the first time will define the start of the observation and the last one
            will be the end of the observation. Note that the higher time resolution
            in `times` the more precise values you will obtain for antenna source
            visibility or determination of the rms noise levels. However, that will
            also imply a longer computing time. Steps of ~5-15 min seem appropriate
            for typical VLBI observations.
        - band : str
            Observing band to conduct the observation. Note that while this is a
            free-format string, it should match the type used when defining the bands
            at which each station can observe, and the default is the str format `XXcm`
            where XX represents the wavelength to observe in cm units.
            You can always check the available bands in `{Stations object}.observing_bands`.
        - datarate : int or astropy.units.Quantity
            Data rate for each antenna. It assumes that all antennas will run at the same
            data rate, which may not be true. If an int is introduce, it will be assumed to
            be given in Mbit/s. Otherwise, you can specify a Quantity with compatible units.
        - subbands : int
            Number of subbands in which the total bandwidth of the observation will be divided
            during correlation.
        - channels : int
            Number of channels for each subband to be created during correlation.
        - polarizations : int (1, 2, 4)
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
        - stations : vlbiplanobs.stations.Stations
            Network of stations that will participate in the given observation.
        - bits : int
            Number of bits at which the data have been recorded (sampled). A typical VLBI observation is
            almost always recorded with 2-bit data.
        """
        self.target = target
        # Because otherwise gstimes is not initialized
        if times is None:
            self._times = None
            self._gstimes = None
        else:
            self.times = times

        self.band = band
        self.datarate = datarate
        self.subbands = subbands
        self.channels = channels
        self.polarizations = polarizations
        self.inttime = inttime
        if stations is not None:
            self.stations = stations
        else:
            self.stations = Stations('empty', [])

        self.bitsampling = bits
        self.ontarget_fraction = ontarget
        self._uv_baseline = None
        self._uv_array = None
        self._rms = None
        self._synth_beam = None


    @property
    def target(self) -> Source:
        """Returns the target source to be observed during the current observation.
        It can return None if the target has not been set yet, showing a warning.
        """
        if self._target is None:
            print("WARNING: 'target' source not set yet but used in 'Observation'.")
        return self._target


    @target.setter
    def target(self, new_target: Source):
        assert isinstance(new_target, Source) or new_target is None, \
                "The new target must be a observation.Source instance, or None."
        self._target = new_target
        # Resets all parameters than depend on the source coordinates
        self._uv_baseline = None
        self._uv_array = None
        self._rms = None
        self._synth_beam = None


    @property
    def times(self) -> Time:
        """Returns the times when the observation runs as an astropy.time.Time object.
        It can return None if the times have not been set yet, showing a warning.
        """
        if self._times is None:
            print("WARNING: 'times' not set yet but used in 'Observation'.")
        return self._times


    @times.setter
    def times(self, new_times):
        assert isinstance(new_times, Time) and len(new_times) >= 2
        self._times = new_times
        self._gstimes = self._times.sidereal_time('mean', 'greenwich')
        self._uv_baseline = None
        self._uv_array = None
        self._rms = None
        self._synth_beam = None


    @property
    def gstimes(self) -> coord.angles.Longitude:
        """Returns the GST times when the observation runs as an astropy.coordinates.angles.Longitude
        object (meaning in hourangle units).
        It can return None if the times have not been set yet, showing a warning.
        """
        if self._gstimes is None:
            print("WARNING: 'times' not set yet but used in 'Observation'.")
        return self._gstimes


    @property
    def duration(self) -> u.Quantity:
        """Returns the total duration of the observation.
        It raises AttributeError if the attribute 'times' has not been set yet.
        """
        if self.times is None:
            raise AttributeError("'times' in 'Observation' has not been set yet.")

        return (self.times[-1] - self.times[0]).to(u.h)


    @property
    def band(self) -> str:
        """Returns the observing band at which the observation will be conducted.
        It can return None if the band has not been set yet, showing a warning.
        """
        if self._band is None:
            print("WARNING: 'band' not set yet but used in 'Observation'.")

        return self._band


    @band.setter
    def band(self, new_band):
        assert (isinstance(new_band, str) and 'cm' in new_band) or (new_band is None)
        self._band = new_band
        self._uv_baseline = None
        self._uv_array = None
        self._rms = None
        self._synth_beam = None


    @property
    def wavelength(self) -> u.Quantity:
        """Returns the central wavelength of the observation.
        """
        assert self.band is not None
        return float(self.band.replace('cm',''))*u.cm


    @property
    def frequency(self) -> u.Quantity:
        """Returns the central frequency of the observations.
        """
        assert self.band is not None
        return 30*u.GHz/self.wavelength.to(u.cm).value


    @property
    def datarate(self) -> u.Quantity:
        """Retuns the data rate (per station) used at which the observation is conducted.
        It can return None if the data rate has not been set yet, showing a warning.
        """
        if self._datarate is None:
            print("WARNING: 'datarate' not set yet but used in 'Observation'.")

        return self._datarate


    @datarate.setter
    def datarate(self, new_datarate):
        """Sets the data rate used at each station during the observation.

        Inputs
        - new_datarate : int or astropy.units.Quantity
            If no units provided, Mbps assumed.
        """
        if new_datarate is None:
            self._datarate = None
        elif isinstance(new_datarate, int):
            self._datarate = new_datarate*u.Mbit/u.s
        elif isinstance(new_datarate, u.Quantity):
            self._datarate = new_datarate.to(u.Mbit/u.s)
        else:
            raise ValueError(f"Unknown type for datarate {new_datarate}" \
                             "(int or astropy.units.Quantity (~bit/s) expected)")

        self._rms = None


    @property
    def subbands(self) -> int:
        """Returns the number of subbands (also known as intermediate frequencies for AIPS users)
        in which the total bandwidth of the observation will be divided during correlation.
        It can return None if the number of subbands has not been set yet, showing a warning.
        """
        if self._subbands is None:
            print("WARNING: 'subbands' not set yet but used in 'Observation'.")

        return self._subbands


    @subbands.setter
    def subbands(self, n_subbands: int):
        assert (isinstance(n_subbands, int) and n_subbands > 0) or n_subbands is None
        self._subbands = n_subbands


    @property
    def channels(self) -> int:
        """Returns the number of channels in which each subband will be divided during correlation.
        It can return None if the number of channels has not been set yet, showing a warning.
        """
        if self._channels is None:
            print("WARNING: 'channels' not set yet but used in 'Observation'.")

        return self._channels


    @channels.setter
    def channels(self, n_channels: int):
        assert (isinstance(n_channels, int) and n_channels > 0) or (n_channels is None)
        self._channels = n_channels


    @property
    def polarizations(self) -> int:
        """Returns the number of polarizations that will be stored in the final data.
        It can return None if the number of polarizations has not been set yet, showing a warning.
        """
        if self._polarizations is None:
            print("WARNING: 'polarizations' not set yet but used in 'Observation'.")

        return self._polarizations


    @polarizations.setter
    def polarizations(self, n_pols: int):
        assert (n_pols in (1, 2, 4)) or (n_pols is None)
        self._polarizations = n_pols


    @property
    def inttime(self) -> u.Quantity:
        """Returns the integration time used when correlating the observation as an astropy.units.Quantity.
        It can return None if the integration time has not been set yet, showing a warning.
        """
        if self._inttime is None:
            print("WARNING: 'inttime' not set yet but used in 'Observation'.")

        return self._inttime


    @inttime.setter
    def inttime(self, new_inttime):
        """Sets the integration time of the observation.
        Inputs
        - new_inttime float/int or astropy.units.Quantity.
            If no units provided, seconds are assumed.
        """
        if new_inttime is None:
            self._inttime = None
        elif isinstance(new_inttime, float) or isinstance(new_inttime, int):
            self._inttime = new_inttime*u.s
        elif isinstance(new_inttime, u.Quantity):
            self._inttime = new_inttime.to(u.s)
        else:
            raise ValueError(f"Unknown type for 'inttime' {new_inttime} (float/int/Quantity(~seconds) expected)")


    @property
    def ontarget_fraction(self) -> float:
        """Fraction of the total observing time spent on the target source.
        It can return None if the ontarget_fraction has not been set yet, showing a warning.
        """
        if self._ontarget is None:
            print("WARNING: 'ontarget_fraction' not set yet but used in 'Observation'.")

        return self._ontarget


    @ontarget_fraction.setter
    def ontarget_fraction(self, ontarget: float):
        assert (0.0 < ontarget <= 1.0) or (ontarget is None)
        self._ontarget = ontarget
        self._rms = None


    @property
    def ontarget_time(self) -> u.Quantity:
        """Total time spent on the target source during the observation.
        It can return None if the ontarget_fraction and duration have not been set yet, showing a warning.
        """
        if (self._ontarget is None):
            raise AttributeError("'ontarget_time' in 'Observation' has not been set yet.")

        return self.duration*self.ontarget_fraction


    @property
    def bandwidth(self) -> u.Quantity:
        """Returns the total bandwidth of the observation.
        It raises AttributeError if the attributes 'polarizations', 'datarate', or 'bitsampling'
        have not been set yet.
        """
        if None in (self.polarizations, self.datarate, self.bitsampling):
            raise AttributeError("'polarizations', 'datarate', and 'bitsampling' in 'Observation'" \
                                 "have not been set yet.")
        pols = self.polarizations % 3 + self.polarizations // 3  # Either 1 or 2
        return (self.datarate/(pols*self.bitsampling*2)).to(u.MHz)


    @property
    def bitsampling(self) -> u.Quantity:
        """Returns the bit sampling at which the data are recorded during the observation.
        """
        return self._bitsampling


    @bitsampling.setter
    def bitsampling(self, new_bitsampling):
        """Sets the bit sampling of the observation.
        Inputs
        - new_bitsampling : int/float or astropy.units.Quantity
            If no units provided, bits are assumed.
        """
        if isinstance(new_bitsampling, float) or isinstance(new_bitsampling, int):
            self._bitsampling = new_bitsampling*u.bit
        elif isinstance(new_inttime, u.Quantity):
            self._bitsampling = new_bitsampling.to(u.bit)
        else:
            raise ValueError(f"Unknown type for {new_bitsampling} (float/int/Quantity(bit) expected)")


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


    def elevations(self) -> dict:
        """Returns the elevation of the target source for each stations participating in the observation
        for all the given observing times.

        Returns
            elevations : dict
                Dictionary where they keys are the station code names, and the values will be
                an astropy.coordinates.angles.Latitude object with the elevation at each observing time.
        """
        elevations = {}
        for a_station in self.stations:
            elevations[a_station.codename] = a_station.elevation(self.times, self.target)
        return elevations


    def altaz(self) -> dict:
        """Returns the altitude/azimuth of the target source for each stations participating
        in the observation for all the given observing times.

        Returns
            altaz : dict
                Dictionary where they keys are the station code names, and the values will be
                an astropy.coordinates.SkyCoord object with the altitude and azimuth
                of the source at each observing time.
        """
        aa = {}
        for a_station in self.stations:
            aa[a_station.codename] = a_station.altaz(self.times, self.target)
        return aa


    def is_visible(self) -> dict:
        """Returns whenever the target source is visible for each station for each time
        of the observation.

        Returns
            is_visible : dict
                Dictionary where they keys are the station code names, and the values will be
                a tuple containing a numpy array with the indexes in the `Observation.times`
                array with the times where the target source can be observed by the station.

                In this sense, you can e.g. call obs.times[obs.is_visible[a_station_codename]]
                to get such times.
        """
        iv = {}
        for a_station in self.stations:
            iv[a_station.codename] = a_station.is_visible(self.times, self.target)
        return iv


    def longest_baseline(self) -> tuple:
        """Returns the longest baseline in the observation.

        Returns
        - ('{ant1}-{ant2}', length) : tuple
            - '{ant1}-{ant2}' : str
                Composed by the codenames of the two antennas (ant1, ant2) conforming the longest baseline.
            - length : astropy.units.Quantity
                The projected length of the baseline as seen from the target source position.
        """
        uv = self.get_uv_baseline()
        longest_bl = {'bl': '', 'value': None}
        for a_bl in uv:
            bl_length = np.sqrt(np.max((uv[a_bl]**2).sum(axis=1)))
            if (longest_bl['value'] is None) or (bl_length > longest_bl['value']):
                longest_bl['bl'] = a_bl
                longest_bl['value'] = bl_length

        return longest_bl['bl'], longest_bl['value']*self.wavelength


    def shortest_baseline(self) -> tuple:
        """Returns the shortest baseline in the observation.

        Returns
        - ('{ant1}-{ant2}', length) : tuple
            - '{ant1}-{ant2}' : str
                Composed by the codenames of the two antennas (ant1, ant2) conforming the shortest baseline.
            - length : astropy.units.Quantity
                The projected length of the baseline as seen from the target source position.
        """
        uv = self.get_uv_baseline()
        shortest_bl = {'bl': '', 'value': None}
        for a_bl in uv:
            bl_length = np.sqrt(np.max((uv[a_bl]**2).sum(axis=1)))
            if (shortest_bl['value'] is None) or (bl_length < shortest_bl['value']):
            # if (bl_length < shortest_bl['value']) or (shortest_bl['value'] is None):
                shortest_bl['bl'] = a_bl
                shortest_bl['value'] = bl_length

        return shortest_bl['bl'], shortest_bl['value']*self.wavelength


    def bandwidth_smearing(self) -> u.Quantity:
        """Returns the bandwidth smearing expected for the given observation.

        The peak response to a point target source decreases at positions farther away from the
        pointing (correlated) sky position due to the frequency averaging performed in the data.

        This function returns the angular separation at which the bandwidth smearing produces
        a reduction of a 10% in the response of the telescope. The field of view should then
        be limited to this range to avoid significant loses.
        """
        return ((49500*u.arcsec*u.MHz*u.km)*self.channels/ \
            (self.longest_baseline()[1]*self.bandwidth/self.subbands)).to(u.arcsec)


    def time_smearing(self) -> u.Quantity:
        """Returns the time smearing expected for the given observation.

        The peak response to a point target source decreases at positions farther away from the
        pointing (correlated) sky position due to the time averaging performed in the data.

        This function returns the angular separation at which the time smearing produces
        a reduction of a 10% in the response of the telescope. The field of view should then
        be limited to this range to avoid significant loses.
        """
        return ((18560*u.arcsec*u.km*u.s/u.cm)* \
                (self.wavelength/(self.longest_baseline()[1]*self.inttime))).to(u.arcsec)


    def datasize(self) -> u.Quantity:
        """Returns the expected size for the output FITS IDI files.

        A regular observation with the European VLBI Network is stored in FITS IDI files,
        typically several 2-GB files. This function provides an estimation of the total
        size for these stored files.
        Note that this function does not take into account down times for the different
        stations. The provided value will thus always be un upper-limit for the real, final,
        value.
        """
        temp = len(self.stations)**2*(self.duration/self.inttime).decompose()
        temp *= self.polarizations*self.subbands*self.channels
        return temp*1.75*u.GB/(131072*3600)


    def thermal_noise(self) -> u.Quantity:
        """Returns the expected rms thermal noise for the given observation.

        Each antenna has a different sensitivity for the observing band (established from
        the SEFD value). This function computes the available baselines at each timestamp
        and estimates the final thermal noise reached when integrating the whole observation.

        Note that this function takes into account when each station can observe the source,
        but does not take into account sensitivity drops doe to external factors like e.g.
        low elevations of the source. The provided thermal noise is also assumed when no
        weighting is applied to the data when imaging. The thermal noise can thus be a bit
        lower if a natural robust weighting is used, but slightly higher if an uniform
        roubst is used.
        """
        # As the determination is computationally heavy, if no parameters have been updated
        # the nit returns the stored value from a previous run.
        if self._rms is not None:
            return self._rms

        main_matrix = np.zeros((len(self.times), len(self.stations)))
        visible = self.is_visible()
        for i,stat in enumerate(self.stations):
            main_matrix[:,i][visible[stat.codename]] = stat.sefd(self.band)
        # Determines the noise level for each time stamp.
        temp = 0.0
        for i,ti in enumerate(main_matrix[:-1]):
            sefds = ti[np.where(ti > 0.0)]
            for j in range(len(sefds)):
                for k in range(j+1, len(sefds)):
                    temp += (self.times[i+1]-self.times[i]).to(u.s).value/(sefds[j]*sefds[k])

        if temp == 0.0:
            # No sources visible
            raise SourceNotVisible('No single baseline can observe the source.')

        temp = 1.0/np.sqrt(temp*self.ontarget_fraction)
        self._rms = ((1.0/0.7)*temp/np.sqrt(self.datarate.to(u.bit/u.s).value/2))*u.Jy
        return self._rms


    def get_uv_baseline(self) -> dict:
        """Returns the (u, v) values for each baseline and each timestamp for which the source
        is visible.

        Returns
        - {'{ant1}-{ant2}': uv_data} : dict
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

        bl_uv_up = {}
        hourangle = (self.gstimes - self.target.ra.to(u.hourangle)).value % 24*u.hourangle
        nstat = len(self.stations)
        # Determines the xyz of all baselines. Time independent
        bl_xyz = np.empty(((nstat*(nstat-1))//2, 3))
        bl_names = []
        s = [ant.location for ant in self.stations]
        for i in range(nstat):
            for j in range(i+1, nstat):
                # An unique number defining a baseline
                k = int( i*(nstat-1) - sum(range(i)) + j-i )
                bl_xyz[k-1,:] = np.array([ii.value for ii in s[i].to_geocentric()]) - \
                             np.array([ii.value for ii in s[j].to_geocentric()])
                bl_names.append("{}-{}".format(self.stations[i].codename,
                                               self.stations[j].codename))

        # Matrix to convert xyz to uvw for each timestamp (but w is not considered)
        m = np.array([[np.sin(hourangle), np.cos(hourangle), np.zeros(len(hourangle))],
                      [-np.sin(self.target.dec)*np.cos(hourangle),
                      np.sin(self.target.dec)*np.sin(hourangle),
                      np.cos(self.target.dec)*np.ones(len(hourangle))]])

        bl_uv = np.array([m[:,:,i] @ bl_xyz.T  for i in range(m.shape[-1])])*u.m
        ants_up = self.is_visible()
        for i,bl_name in enumerate(bl_names):
            ant1, ant2 = bl_name.split('-')
            bl_up = (np.array([a for a in ants_up[ant1][0] if a in ants_up[ant2][0]]), )
            if len(bl_up[0]) > 0:
                bl_uv_up[bl_name] = (bl_uv[:,:,i][bl_up]/self.wavelength).decompose()

        if len(bl_uv_up.keys()) == 0:
            raise SourceNotVisible

        self._uv_baseline = bl_uv_up
        return bl_uv_up


    def get_uv_array(self) -> np.ndarray:
        """Returns the (u, v) values for each baseline and each timestamp for which the source
        is visible.

        The difference with `get_uv_baseline` is that `get_uv_array` only returns the (u,v)
        values, dropping the information of baselines and times to which these values belong to.

        Returns a (N, 2)-dimensional numpy.ndarray containing all N (u,v) points resulting for
        each timestamp and baseline. The (u,v) values are given in lambda units.
        Note that complex conjugate values are not provided.

        Exceptions
        - It may raise the exception SourceNotVisible if no baselines can observe the source
          at all during the observation.
        """
        if self._uv_array is not None:
            return self._uv_array

        bl_uv_up = self.get_uv_baseline()
        tot_length = 0
        for bl_name in bl_uv_up:
            tot_length += bl_uv_up[bl_name].shape[0]

        uvvis = np.empty((tot_length, 2))
        i = 0
        for bl_name in bl_uv_up:
            uvvis[i:i+bl_uv_up[bl_name].shape[0],:] = bl_uv_up[bl_name]
            i += bl_uv_up[bl_name].shape[0]
        self._uv_array = uvvis
        return self._uv_array


    def synthesized_beam(self) -> dict:
        """Estimates the resulting synthesized beam of the observations based on
        the expected (u,v) coverage.

        This is just an estimation made by a ellipse fitting to the (u, v) coverage,
        from which we obtain the resolution on the two axes following
        https://science.nrao.edu/facilities/vlba/docs/manuals/oss/ang-res
            theta_HPBW (mas)  \sim 2063 x lambda(cm)/b_max^km

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

        resolution = lambda bl : ((2.063e8*u.mas)/bl).to(u.mas)
        uvvis = self.get_uv_array()
        # Transform the uv points into r,theta (polar) points
        uvvis_polar = np.empty_like(uvvis)
        uvvis_polar[:,0] = np.sqrt((uvvis**2).sum(axis=1)) # radius
        uvvis_polar[:,1] = np.arctan2(uvvis[:,1], uvvis[:,0]) # theta
        # Defines the BMAJ and PA
        bl_bmaj = np.max(uvvis_polar[:,0])
        bl_bmaj_theta = uvvis_polar[:,1][np.where(uvvis_polar[:,0] == bl_bmaj)][0]
        # Gets the BMIN and an orthogonal projection
        bl_bmin_theta = ( bl_bmaj_theta + np.pi/2 ) % (2*np.pi)
        bl_bmin = np.max(np.abs(uvvis.dot(np.array([np.cos(bl_bmin_theta),
                                                    np.sin(bl_bmin_theta)]))))

        self._synth_beam = {'bmaj': resolution(bl_bmin), 'bmin': resolution(bl_bmaj),
                            'pa': (bl_bmaj_theta*u.rad).to(u.deg)}
        return self._synth_beam

    # def get_dirtymap(self):
    #     uvdata = self._get_uv()
    #     # Generates a N 2-D array with all uv data.
    #     tot_length = 0
    #     for bl_name in uvdata:
    #         tot_length += uvdata[bl_name].shape[0]
    #
    #     uvvis = np.empty((tot_length, 2))
    #     i = 0
    #     for bl_name in uvdata:
    #         uvvis[i:uvdata[bl_name].shape[0]] = uvdata[bl_name]
    #         i += uvdata[bl_name].shape[0]


    def print_obs_times(self, date_format='%d %B %Y'):
        """Returns the time range (starttime-to-endtime) of the observation in a smart way.
        If the observation lasts for less than one day it omits the end date:
                20 January 1971 10:00-20:00 UTC
                GST range: 05:00-15:00

        If the observation ends the day after, then it returns:
                20 January 1971 10:00-20:00 UTC (+1d)
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






