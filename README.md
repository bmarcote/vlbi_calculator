> Thank you everyone for the wait!  PlanObs is back with tons of new features.
> That also imply that for the time being it may contain some minor bugs, and a few features may be missing. They will be back incrementally during the upcoming updates.


# EVN Observation Planner


The EVN Observation Planner is a tool to determine the visibility of a given astronomical source when planning very-long-baseline-interferometry (VLBI) observations. The tool is specially written for the preparation of observations with the [European VLBI Network (EVN)](https://www.evlbi.org), but it can be used for any kind of VLBI observations than can be currently arranged (e.g. with the [Very Long Baseline Array, VLBA](https://public.nrao.edu/telescopes/vlba/); [the Australian Long Baseline Array, LBA](https://www.atnf.csiro.au/vlbi/overview/index.html); [eMERLIN](http://www.merlin.ac.uk/e-merlin/index.html); or [the global mm-VLBI array](https://www3.mpifr-bonn.mpg.de/div/vlbi/globalmm/), for example). An ad-doc VLBI array can also be quickly configured.

In addition to the determination of the source visibility by the different antennas, the EVN Observation Planner would provide an estimation of the expected rms noise level (sensitivity) reached during the planned observations, and an estimation of the resolution. The EVN Observation Planner can thus be used while [preparing an observing proposal](https://www.evlbi.org/using-evn).
Note that the EVN Observation Planner has been designed as a more complete version of the previous [EVN Calculator](http://old.evlbi.org/cgi-bin/EVNcalc.pl).

**Documentation:** <https://bmarcote.github.io/vlbi_calculator/> (installation, all command-line modes, scheduling, and the Python API). The list of changes per version is in [CHANGES.txt](CHANGES.txt).



## It runs both online and locally!

You can make use of the EVN Observation Planner just by going to [the online tool hosted at JIVE](https://planobs.jive.eu), without installing anything.


But if you want to run it in your local machine, you can also install the package via `pip` (Python 3.12 or newer is required):

```bash
python3 -m pip install vlbiplanobs
```


Once you have it installed, you can simply run it by typing `planobs server` in the terminal.  It will start to run the server and you will be able to access it in your browser by following the typed url (by default http://127.0.0.1:8050/).


> **But PlanObs also has a lovely command-line interface!**

The EVN Observation Planner can also be used through the terminal or inside your own Python program without the need of running a server.


### Command-line interface (CLI)

_This is great for when you want to plan, or verify the feasibility of, a VLBI observation quickly, or you have multiple sources to observe already defined._


You only need to type:

```bash
planobs
```


In your terminal to encounter all options.


**Example cases**

1. You want to know when a source (e.g. 'Altair') can be observed by a given VLBI network (e.g. the EVN and eMERLIN). It is as easy as:

```bash
planobs -b 6cm -t 'Altair' --network EVN eMERLIN
```


If you want to specify an epoch and some particular stations, you can do:

```bash
planobs --band 6cm --target 'Altair' --stations Ef Hh Ir Mc Tr T6 O8 Wb Cm --epoch '2020-06-15 20:00' --duration 12
```

The epoch is the start of the observation in UTC, and the duration is given in hours (decimals are fine, e.g. `-d 1.5`).
A target can be a source name, coordinates (`'hh:mm:ss dd:mm:ss'`), or `'name/coordinates'`, and several targets can
be given at once. Use `--list-networks`, `--list-antennas` and `--list-bands` to see the available values.

Add `-o report.pdf` (or `.txt`, `.md`, `.json`) to save all inputs and results to a file, and `--sched EXPCODE` to
produce a SCHED `.key` schedule file (with `--phasecal`, `--check-source`, `--fringefinders` and `--polcal` to include
calibrators, and `--template`/`--get-key-template` to use your own `.key` template).

2. Other modes are available as subcommands (run `planobs <mode> -h` for details):

```bash
planobs fringefinders -n EVN -e '2025-06-15 20:00' -d 8 -b 6cm   # bright fringe finders
planobs phasecals 'M87' -b 6cm                                  # phase calibrators near a target
planobs source '3C273'                                          # source information
planobs antenna Ef                                              # antenna information
planobs server                                                  # local web GUI
```

The same option names are used in every mode: `-b/--band`, `-n/--network`, `-s/--stations`, `-e/--epoch`,
`-d/--duration`, `-t/--target`, and `-l/--max-lines`. The option names used before v5.1.0 (`-t1`, `--starttime`,
`--targets`, `--n-sources`, `--catalog-file`) still work but print a deprecation warning.

The full documentation is [available online](https://bmarcote.github.io/vlbi_calculator/), and its sources are in the
`docs/` directory.



## Station additions

The information about each station is stored in an independent file under `src/vlbiplanobs/data/stations_catalog.inp` (following a [Python configuration file format](https://docs.python.org/3/library/configparser.html)). Then, the addition or update of a new station is extremely straightforward. You can manually add a new station by introducing a new entry in the file with the following fields and syntax:

```ini
[Station Name]

station = Station Name
code = # An unique code to identify the station (typically two letters).
networks = # A comma-separated list of the VLBI networks that the station can join to observe (it can be empty).
country =  # Country where the station is located.
diameter =  # Station diameter in free format (e.g. '30 m', or '30 x 20 m' is often used for interferometers composed of 30 20-m antennas).
position = X, Y, Z  # Geocentric coordinates of the station, in meters.
real_time = no/yes  # If the station can participate in real-time correlation observations (e.g. e-EVN).
sefd_YY = ZZ   # Estimated System Equivalent Flux Density (SEFD) of the station (ZZ, in Jy) at the observing wavelength YY in cm. There should be one entry per observing band.
# The lack of an entry means that the station is unable to observe at such band.

# Optional fields:
sched_name = # Name of the station in SCHED/pySCHED (by default, the name of the entry).
mount = ALTAZ  # Mount type (ALTAZ or EQUAT). Together with ax1rate, ax2rate (slewing speeds, deg/s), ax1lim, ax2lim
# (axis limits 'min, max', deg), and optionally ax1acc, ax2acc (deg/s/s), it is used to compute the slewing times.
horizon_az = 0, 90, 180, 270, 360  # Azimuths (deg) of the local horizon mask...
horizon_el = 10, 12, 10, 15, 10    # ...and the minimum elevation (deg) at each of them (same number of values).
maxdatarate = 2048  # Maximum data rate (Mbps), if lower than the one of the network.
maxdatarate_512 = EVN  # Or a maximum data rate (here 512 Mbps) that only applies within the given network(s).
group = # Label to group several entries in the GUI (e.g. the different configurations of an array).
decommissioned = yes  # If the station is no longer operational.
```

Note that the comments after each value are only explanatory here: the catalog does not support comments in the same line as a value (only in separate lines starting with `#`).

You can also keep your own catalog in a separate file and use it with `planobs --station-catalog FILE`.


We are more than glad to integrate any additional station that can be relevant to the purposes of this program.

If you have any suggestion, please open an issue in the GitHub repository, or [contact the author](mailto:marcote@jive.eu).
