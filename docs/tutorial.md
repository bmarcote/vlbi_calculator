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


## Finding Fringe Finder Sources

For VLBI observations, you need bright calibrator sources to find and correct instrumental delays. These are called fringe finder sources. You can find suitable ones using:

```bash
planobs_fringefinder -s Ef Hh Mc Tr -t '2025-03-15 08:00' -d 8 -b 6cm
```

This command will:
- Search the RFC catalog for bright calibrator sources
- Check visibility for the specified stations during the observation time
- Show flux densities, elevation information, and links to detailed source data

### Fringe Finder Options

- `-s/--stations`: List of antenna codenames that will participate
- `-t/--starttime`: Start time in 'YYYY-MM-DD HH:MM' format (UTC)
- `-d/--duration`: Duration of the observation in hours
- `-b/--band`: Observing band (e.g., '18cm', '21cm', '13cm', '6cm', '5cm', '3.6cm', '2cm', '1.3cm', '0.7cm')
- `--min-flux`: Minimum unresolved flux threshold in Jy (default: 0.5)
- `--min-elevation`: Minimum elevation in degrees (default: 20)
- `-l/--max-lines`: Maximum number of sources to return (default: 20)
- `--require-all`: Require source to be visible by ALL stations (default: False)

## Finding Phase Calibrators

Phase referencing requires finding bright calibrator sources near your target source. You can find these using:

```bash
planobs_phasecal -t 'M87' -b 6cm --max-separation 5 --min-flux 0.1
```

This command will:
- Search for calibrator sources near your target
- Show angular separation from the target
- Display flux information at the specified band
- Provide links to detailed source information

### Phase Calibrator Options

- `-t/--target`: Target source name (J2000 or IVS name from RFC catalog)
- `-b/--band`: Observing band (optional, checks all bands if not provided)
- `--max-separation`: Maximum angular separation in degrees (default: 5.0)
- `--min-flux`: Minimum unresolved flux threshold in Jy (default: 0.1)
- `-n/--n-sources`: Maximum number of sources to return (default: all)
- `--catalog-file`: Path to custom RFC catalog file

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

