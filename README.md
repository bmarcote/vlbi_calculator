# EVN Observation Planner


The EVN Observation Planner is a tool to determine the visibility of a given astronomical source when planning very-long-baseline-interferometry (VLBI) observations. The tool is specially written for the preparation of observations with the [European VLBI Network (EVN)](https://www.evlbi.org), but it can be used for any kind of VLBI observations than can be currently arranged (e.g. with the [Very Long Baseline Array, VLBA](https://public.nrao.edu/telescopes/vlba/); [the Australian Long Baseline Array, LBA](https://www.atnf.csiro.au/vlbi/overview/index.html); [eMERLIN](http://www.merlin.ac.uk/e-merlin/index.html); or [the global mm-VLBI array](https://www3.mpifr-bonn.mpg.de/div/vlbi/globalmm/), for example). An ad-doc VLBI array can also be quickly configured.

In addition to the determination of the source visibility by the different antennas, the EVN Observation Planner would provide an estimation of the expected rms noise level (sensitivity) reached during the planned observations, and an estimation of the resolution. The EVN Observation Planner can thus be used while [preparing an observing proposal](https://www.evlbi.org/using-evn).
Note that the EVN Observation Planner has been designed as a complementary, but more featured, version of the current [EVN Calculator](http://old.evlbi.org/cgi-bin/EVNcalc.pl).



## It runs online!

You can make use of the EVN Observation Planner just by going to [the Heroku running server](https://vlbi-calculator.herokuapp.com/), without installing anything. Note that this is a temporary link, and soon it will be hosted in its final location inside the [EVN](https://www.evlbi.org) website.


It only requires the minimal information to be able to compute the results of the observation:

- Select the observing band (at which frequency, or wavelength, do you want to observe?).
- Select a default VLBI network, or make an ad-hoc one by selecting manually the worldwide antennas that you want to use. You will see that only the antennas that can observe at the given band would be selectable.
- Enter the start and end of the planned observation, and the coordinates (in J2000 format) of the main source you want to observe.
- You can still tune more technical details of the observation, like the expected data rate, number of subbands, channels, or integration time. But default values will always been set automatically for your help.

After all this, you only need to press a button (`computer observation`) and you will presented with a detailed output in the different tabs:

- A summary showing the size of the data that you may expect to be downloaded (once correlated), the longest and shortest baselines in the array, the expected synthesized beam (resolution), and rms noise level of the resulting image, and the limitations in the field of view due to the frequency and time averaging.
- A couple of plots showing the elevation of the source for all the selected antennas, and the time ranges when they can observe the source.
- A plot showing the expected _u, v_ coverage. Note that depending on how filled the (_u,v_) plane is, the better reconstructed the resulting image will be.


You can also run it locally by downloading the code and running:

```bash
python app.py
```



## But you can also use it inside your Python programs or interactively!

The EVN Observation Planner can also be used inside a Python environment or inside your own programs without the need of running a server.

*Note that this part is still under development and a detailed use will be specified soon*





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




