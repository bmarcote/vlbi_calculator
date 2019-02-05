#util functions for seffers
import numpy as np
import datetime as dt

from astropy import coordinates as coord
from astropy.time import Time
from astropy import units as u


# useful functions

def get_coordinates(text_coord):
    """Given a string representing coordinates in the form hh:mm:ss dd:mm:ss
    if returns the astropy.coordinates from this object.
    """
    ra, dec = text_coord.strip().split(' ')
    ra = [float(i) for i in ra.split(':')]
    dec = [float(i) for i in dec.split(':')]
    return coord.SkyCoord('{:n}h{:n}m{}s {:n}d{:n}m{}s'.format(*ra, *dec))


def get_time(text_time):
    """Returns a datetime.datetime object representing the input text_time that
    should have the form DD/MM/YYYY HH:MM
    """
    # return Observer.datetime_to_astropy_time(dt.datetime.strptime(text_time, '%d/%m/%Y %H:%M'))
    the_time = dt.datetime.strptime(text_time, '%d/%m/%Y %H:%M')
    return Time(the_time.strftime('%Y-%m-%d %H:%M'))
    #date = [int(i) for i in date.split('/')]


def get_obs_times(start_time, duration, interval=0.2):
    """Returns the times of the observation from the starting to the end
    Duration and interval must be in hours, without units.
    """
    return start_time + np.arange(0.0, duration+interval/2., interval)*u.h


def times_elev_with_source_above(source_coordinates, antennas, obs_times):
    pass



def get_thermal_noise(antennas, times, data_rate, band):
    pass


def get_bandwidth_smearing():
    pass


