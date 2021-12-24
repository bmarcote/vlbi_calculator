"""This is a script that loads a dummy observation
(using all parts without the graphical, web-style, interface).

By running it in an interactive mode it allows me to test all parts.
Note that this is not a full test...
"""
import os
import sys
import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.dates as mdates
from astropy import coordinates as coord
from astropy import units as u
from astropy.time import Time
from astropy.io import ascii
from vlbiplanobs import freqsetups as fs
from vlbiplanobs import stations
from vlbiplanobs import observation




target = observation.Source('10h58m29.6s +81d33m58.8s', name='Target')

obs = observation.Observation(target=target)
obs.times = Time('2020-06-15 20:00', scale='utc') + np.arange(0, 720, 10)*u.min
obs.band = '18cm'
obs.datarate = 1024
obs.subbands = 8
obs.channels = 32
obs.polarizations = 2
obs.inttime = 2

all_stations = stations.Stations.get_stations_from_configfile()

def get_selected_antennas(list_of_selected_antennas):
    """Given a list of antenna codenames, it returns a Stations object containing
    all given antennas.
    """
    selected_antennas = stations.Stations('Observation', [])
    for ant in list_of_selected_antennas:
        selected_antennas.add(all_stations[ant])
    return selected_antennas



# evn6 = ['Ef', 'Jb2', 'On', 'Hh', 'T6', 'Wb', 'Sv', 'Zc']
evn6 = ['Ef', 'Jb2', 'On', 'Hh', 'T6', 'Wb', 'Sv', 'Zc', 'Pa', 'Mp', 'Ho', 'Nl', 'Pt', 'Sc', 'Kp', 'Hn']
obs.stations = get_selected_antennas(evn6)


elevs = obs.elevations()
srcup = obs.is_visible()
beam = obs.synthesized_beam()
rms = obs.thermal_noise()
uvdata = obs.get_uv_array()
dirty_map_natural = obs.get_dirtymap(robust='natural')
dirty_map_uniform = obs.get_dirtymap(robust='uniform')
# fig, ax = plt.subplots()


# for ant in elevs:
#     ax.plot(obs.times.datetime[srcup[ant]], elevs[ant][srcup[ant]].deg, '-', label=ant)


# ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
# ax.set_xlabel(f"Time (UTC) - {obs.times.datetime[0].strftime('%d-%m-%Y')}")







def get_decimal_hour(dts):
    return np.array([t.hour+(t.minute/60)+(t.second/3600) for t in dts])


def t2gst(t):
    # Expects a float from the input
    # t2 = np.array([tt.hour+tt.minute/60+tt.second/3600 for tt in t])
    gst_t_slope = 1.00273791 # Is this universal? it is about 16.4 min
    # get the offset from the first timestamp
    offset = -obs.times.datetime[0].hour-obs.times.datetime[0].minute/60 + \
             obs.gstimes[0].hour
    return (t*gst_t_slope + offset) % 24
    # TODO: I see a 0.027379 hour offset in the example... is this stable?

def gst2t(t):
    gst_t_slope = 1.00273791 # Is this universal? it is about 16.4 min
    # get the offset from the first timestamp
    offset = -obs.times.datetime[0].hour-obs.times.datetime[0].minute/60 + \
             obs.gstimes[0].hour
    return ((t - offset)/gst_t_slope) % 24



# ax_gst = ax.secondary_xaxis("top", functions=(t2gst, gst2t))
# ax_gst.cla()






