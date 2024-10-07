**WARNING**:
_I would like to note that I am currently developing a major new version (v2.0) of the Planner. Better UI, much faster, and with a better API. Meanwhile, most of the code is broken in the current repo. Please download the last stable version if you want to use it. But if you have some wishes, I will be happy to receive feedback!_

# EVN Observation Planner


The EVN Observation Planner is a tool to determine the visibility of a given astronomical source when planning very-long-baseline-interferometry (VLBI) observations. The tool is specially written for the preparation of observations with the [European VLBI Network (EVN)](https://www.evlbi.org), but it can be used for any kind of VLBI observations than can be currently arranged (e.g. with the [Very Long Baseline Array, VLBA](https://public.nrao.edu/telescopes/vlba/); [the Australian Long Baseline Array, LBA](https://www.atnf.csiro.au/vlbi/overview/index.html); [eMERLIN](http://www.merlin.ac.uk/e-merlin/index.html); or [the global mm-VLBI array](https://www3.mpifr-bonn.mpg.de/div/vlbi/globalmm/), for example). An ad-doc VLBI array can also be quickly configured.

In addition to the determination of the source visibility by the different antennas, the EVN Observation Planner would provide an estimation of the expected rms noise level (sensitivity) reached during the planned observations, and an estimation of the resolution. The EVN Observation Planner can thus be used while [preparing an observing proposal](https://www.evlbi.org/using-evn).
Note that the EVN Observation Planner has been designed as a complementary, but more featured, version of the current [EVN Calculator](http://old.evlbi.org/cgi-bin/EVNcalc.pl).



## It runs online!

You can make use of the EVN Observation Planner just by going to [the online tool hosted at JIVE](https://planobs.jive.eu), without installing anything.


It only requires the minimal information to be able to compute the results of the observation:

- Select the observing band (at which frequency, or wavelength, do you want to observe?).
- Select a default VLBI network, or make an ad-hoc one by selecting manually the worldwide antennas that you want to use. You will see that only the antennas that can observe at the given band would be selectable.
- Enter the start and end of the planned observation, and the coordinates (in J2000 format) of the main source you want to observe.
- You can still tune more technical details of the observation, like the expected data rate, number of subbands, channels, or integration time. But default values will always been set automatically for your help.

After all this, you only need to press a button (`computer observation`) and you will presented with a detailed output in the different tabs:

- A summary showing the size of the data that you may expect to be downloaded (once correlated), the longest and shortest baselines in the array, the expected synthesized beam (resolution), and rms noise level of the resulting image, and the limitations in the field of view due to the frequency and time averaging.
- A couple of plots showing the elevation of the source for all the selected antennas, and the time ranges when they can observe the source.
- A plot showing the expected _u, v_ coverage. Note that depending on how filled the (_u,v_) plane is, the better reconstructed the resulting image will be.



## Installing it locally

But if you want to run it in your local machine, you can also install the package simply by running

```bash
python setup.py install
```

Soon we will upload it to PyPy so you will be able to install it from `pip`.


Note that the current version requires the package `astropy` **version 4.0.1** and latest `astroplan`. This restriction in the version of `astropy` is produced by a bug in versions <4.0.1 only triggered when multiple instances run the program at the same time [(see issue #10114 from astropy)](https://github.com/astropy/astropy/issues/10114). If you are running vlbiplanobs only through the command line you will not be affected. On the other hand, versions >4.0.1 are not supported by the current version of `astroplan` (0.6). This will be fixed once version 0.7 is released, unlocking the more recent versions of `astropy`.


Once you have it installed, you can simply run it by typing `vlbiplanobs` in the terminal.  It will start to run the server and you will be able to access it in your browser by following the typed url (likely http://0.0.0.0:8050/).




## But you can also use it inside your Python programs or interactively!

The EVN Observation Planner can also be used inside a Python environment or inside your own programs without the need of running a server, ignoring the Dash server.


```python
# Packages used in the example:
import numpy as np
from astropy import units as u

import vlbiplanobs

# Two main modules that can also be imported directly
from vlbiplanobs import stations
from vlbiplanobs import observation

# You can define your source and telescopes before setting the observation:

# You can define the source you want to observe with:
source = observation.Source(coordinates='XXhXXmXXs XXdXXmXXs', name='my_source')

# To can then retrieve the coordinates with:
# (as in a astropy.coordinates.angles.Longitude/Latitude object)
source.ra
source.dec
# Or retrieve the full astropy.coordinates object as
source.coord


# You can import all stations that are known by default:
all_stations = stations.Stations.get_stations_from_configfile()
# Note that you can get the codenames with
all_stations.keys()

# Or just a subset of them by selecting manually
# the stations you want from their codenames (e.g.):
my_stations = stations.Stations.get_stations_from_configfile(codenames=('Ef', 'Ys', 'Wb'))

# Where you can also select the antennas that can observe at a given band:
stations18cm = my_stations.stations_with_band('18cm')


# Finally, you can set the observation
obs = observation.Observation(target=source,
    times=observation.Time('1967-04-17 10:00') + np.arange(0, 600, 15)*u.min,  # list of times covering the observation.
    band='18cm',  # must be a string with the format XXcm.
    datarate=1024, # Mbps
    subbands=8, # no. subbands (or IFs).
    channels=64, # no. of channels per subband
    polarizations=4, # no. of polarizations (1: single, 2: dual, 4: full polarization)
    inttime=2,  # integration time in seconds (or astropy.Quantity)
    ontarget=0.7, # fraction of the total observing time spent in the target source (affects to the estimated noise level)
    stations=stations18cm,   # add a Stations object containing all stations that will observe
    bits=2)  # no. of bits used for data recording (2 in typical VLBI observations)

# You can get the wavelength or frequency of the observation
print(obs.wavelength)
print(obs.frequency)

# The total bandwidth of the observation
print(obs.bandwidth)
# Or per subband:
print(obs.bandwidth/obs.subbands)

# Get the times (UTC) in datetime format, or GST:
print(obs.times.datetime)
print(obs.gstimes)


# Obtain the source elevation along the observation for each antenna
obs.elevations()  # Returns a dict with the codename of the station as key and an numpy.array with the elevations as value.

# or the altitude/azimuth of the source, following the same approach
obs.altaz()

# If you just want to know at which times the source will be visible per station:
obs.is_visible()

# Obtain the expected thermal rms noise for the whole observation
obs.thermal_noise()

# And the expected synthesized beam (using a neutral weighting of zero)
obs.synthesized_beam()  # returns a dict with 'bmaj', 'bmin', 'pa'.

# And some additional useful information as (but not limited to):
obs.longest_baseline()  # returning ((antenna1,antenna2), baseline_length)
obs.shortest_baseline()
obs.bandwidth_smearing()
obs.time_smearing()

```

You can then use of favourite tools or your own scripts to process this information.



## Station additions

The information about each station is stored in an independent file under `data/stations_catalog.inp` (following a [Python configuration file format](https://docs.python.org/3/library/configparser.html)). Then, the addition or update of a new station is extremely straightforward. You can manually add a new station by introducing a new entry in the file with the following fields and syntax:

```python
[Station Name]

station = Station Name
code = # An unique code to identify the station (it can be the abbreviation of the full station name).
network = # If the station belongs to one of the known VLBI Networks, or 'Other' otherwise.
possible_networks = # a comma-separated list of possible VLBI Networks that the station can join to observe.
country =  # Country where the station is located.
diameter =  # station diameter in free format (e.g. '30 m' or '30 x 20 m' is often used for the case of interferometers composed of 30 20-m antennas).
position = X, Y, Z  # Geocentric coordinates of the station.
min_elevation = XX  # minimum elevation the station can observe, in degrees. By default it is 10 deg if not specified.
real_time = no/yes  #  In case the station can participate in real-time correlation observations (e.g. e-EVN). By default 'no'.
SEFD_YY = ZZ   # Multiple inputs providing the estimated System Equivalent Flux Density (SEFD) of the station (ZZ measured in Jy) at the observing wavelength YY in cm. There should be one entry per observing band.
# The lack of entries would be understood as the station is unable to observe at such band.
```


We are more than glad to integrate any additional station that can be relevant to the purposes of this program.

If you have any suggestion, please open an issue in the GitHub repository, or [contact the author](mailto:marcote@jive.eu).






