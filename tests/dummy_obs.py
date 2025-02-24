"""This is a script that loads a dummy observation
(using all parts without the graphical, web-style, interface).

By running it in an interactive mode it allows me to test all parts.
Note that this is not a full test...
"""
import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.dates as mdates
from astropy import units as u
from astropy.time import Time
from vlbiplanobs import stations as stats
from vlbiplanobs import observation as obs
from vlbiplanobs import sources as src

print('Getting all stations')
all_stations = stats.Stations()

print('Selecting EVNs')
print(f"With only defaults: {all_stations.filter_networks('EVN', only_defaults=True)}")
print(f"With all: {all_stations.filter_networks('EVN', only_defaults=False)}")
print(f"With eMERLINs (defaults only): {all_stations.filter_networks(['EVN', 'eMERLIN'], only_defaults=True)}")
print(f"As string?: {all_stations.filter_networks('EVN,eMERLIN', only_defaults=True)}")

target = src.Source('Target', '10h58m29.6s +81d33m58.8s')
pcal = src.Source('PCAL', '10h59m29.6s +81d33m28.8s')
ccal = src.Source('CHECK', '10h58m00.6s +81d34m28.8s')
scans = src.ScanBlock([src.Scan(source=target, duration=3.5*u.min), src.Scan(source=pcal, duration=1.5*u.min),
                       src.Scan(source=ccal, every=3)])

obs = obs.Observation(band='6cm', stations=all_stations.filter_networks(['EVN', 'eMERLIN'], only_defaults=True),
                      scans=scans, times=Time('2020-06-15 20:00', scale='utc') + np.arange(0, 720, 10)*u.min,
                      datarate=1024*u.Mbit/u.s, subbands=8, channels=64, polarizations=4, inttime=2*u.s)


evn6 = ['Ef', 'Jb2', 'On', 'Hh', 'T6', 'Wb', 'Sv', 'Zc']

# elevs = obs.elevations()
# srcup = obs.is_visible()
# beam = obs.synthesized_beam()
# rms = obs.thermal_noise()
# uvdata = obs.get_uv_array()
# dirty_map_natural = obs.get_dirtymap(robust='natural')
# dirty_map_uniform = obs.get_dirtymap(robust='uniform')
# # fig, ax = plt.subplots()


# for ant in elevs:
#     ax.plot(obs.times.datetime[srcup[ant]], elevs[ant][srcup[ant]].deg, '-', label=ant)


# ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
# ax.set_xlabel(f"Time (UTC) - {obs.times.datetime[0].strftime('%d-%m-%Y')}")


def get_decimal_hour(dts):
    return np.array([t.hour+(t.minute/60)+(t.second/3600) for t in dts])


def t2gst(t):
    # Expects a float from the input
    # t2 = np.array([tt.hour+tt.minute/60+tt.second/3600 for tt in t])
    gst_t_slope = 1.00273791  # Is this universal? it is about 16.4 min
    # get the offset from the first timestamp
    offset = -obs.times.datetime[0].hour-obs.times.datetime[0].minute/60 + obs.gstimes[0].hour
    return (t*gst_t_slope + offset) % 24
    # TODO: I see a 0.027379 hour offset in the example... is this stable?


def gst2t(t):
    gst_t_slope = 1.00273791  # Is this universal? it is about 16.4 min
    # get the offset from the first timestamp
    offset = -obs.times.datetime[0].hour - obs.times.datetime[0].minute/60 + obs.gstimes[0].hour
    return ((t - offset)/gst_t_slope) % 24


# ax_gst = ax.secondary_xaxis("top", functions=(t2gst, gst2t))
# ax_gst.cla()
