# Tutorial



## Quick run

1. You want to know when a source (e.g. 'Altair') can be observed by a given VLBI network (e.g. the EVN and eMERLIN). It is as easy as:

```bash
planobs -b 6cm -t 'Altair' --network EVN eMERLIN
```


If you want to specify an epoch and some particular stations, you can do:

```bash
planobs --band 6cm --target 'Altair' --stations Ef Hh Ir Mc Tr Hh T6 O8 Wb Cm --starttime '2020-06-15 20:00' --duration 12
    --no-gui   # if you don't want to visualize the plots in the browser but only within the terminal
```












## Customizing the stations' catalog


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

