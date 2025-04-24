> Thank you everyone for the wait!  PlanObs is back with tons of new features.
> That also imply that for the time being it may contain some minor bugs, and a few features may be missing. They will be back incrementally during the upcoming updates.


# EVN Observation Planner


The EVN Observation Planner is a tool to determine the visibility of a given astronomical source when planning very-long-baseline-interferometry (VLBI) observations. The tool is specially written for the preparation of observations with the [European VLBI Network (EVN)](https://www.evlbi.org), but it can be used for any kind of VLBI observations than can be currently arranged (e.g. with the [Very Long Baseline Array, VLBA](https://public.nrao.edu/telescopes/vlba/); [the Australian Long Baseline Array, LBA](https://www.atnf.csiro.au/vlbi/overview/index.html); [eMERLIN](http://www.merlin.ac.uk/e-merlin/index.html); or [the global mm-VLBI array](https://www3.mpifr-bonn.mpg.de/div/vlbi/globalmm/), for example). An ad-doc VLBI array can also be quickly configured.

In addition to the determination of the source visibility by the different antennas, the EVN Observation Planner would provide an estimation of the expected rms noise level (sensitivity) reached during the planned observations, and an estimation of the resolution. The EVN Observation Planner can thus be used while [preparing an observing proposal](https://www.evlbi.org/using-evn).
Note that the EVN Observation Planner has been designed as a more complete version of the previous [EVN Calculator](http://old.evlbi.org/cgi-bin/EVNcalc.pl).



## It runs both online and locally!

You can make use of the EVN Observation Planner just by going to [the online tool hosted at JIVE](https://planobs.jive.eu), without installing anything.


But if you want to run it in your local machine, you can also install the package via `pip`:

```bash
python3 -m pip install vlbiplanobs
```


Once you have it installed, you can simply run it by typing `planobs --server` in the terminal.  It will start to run the server and you will be able to access it in your browser by following the typed url (likely http://0.0.0.0:8050/).


> **But PlanObs also has a lovely command-line interface!

The EVN Observation Planner can also be used through the terminal or inside your own Python program without the need of running a server.


### Command-line interface (CLI)

_This is great for when you want to plan, or verify the feasibility of, a VLBI observation quickly, or you have multiple sources to observe already defined._


You only need to type:

```bash
planobs
```


In your terminal to encounter all options.


**Example cases**

1. You want to know when a source (e.g. 'Altair') can be observed by a given VLBI network (e.g. the EVN). It is as easy as:

```bash
planobs --network EVN --target Altair ........
    --no-gui   # if you don't want to visualize the plots in the browser but only within the terminal
```



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






