
# EVN Source Elevation App


## Program structure

- app.py
- README.md
- docs/
- src/
    - stations.py
    - sources.py
    - plots.py
    - functions.py
- data/
    - antenna_positions.txt
    - antenna_sefd.txt
    - frequency_setups.py # or .ini. containting bands, data_rates, etc


### app.py

Layout:
    - band  (keyval)
    - array
    - starttime
    - endtime
    - source
    - onsourcetime
    - datarate
    - subbands
    - channels
    - pols
    - inttime


### stations.py

Station
    - observer : Observer
    - name : str
    - codename : str
    - network : str
    - location : EarthLocation
    - bands : list
    - sefds : dict
    - min_elevation
    + __init__(name, codename, location, freqs_sefds, min_elevation=20)
    + elevation(obs_times, target) --> ndarray
    + is_visible(obs_times, target) --> ndarray of bool?
    + has_band(band) --> bool
    + sefd(band) --> float

Stations
    - name
    - stations
    - number_of_stations
    + __init__(name, stations)
    + add(a_station)
    + keys()  : tuple

    + stations_with_freq()



### observation.py

Source
    - name : str
    - coord : coord.SkyCord
    - dec : coord.Latitute
    - ra : coord.Longitude
    + __init__(coordinates : coord.SkyCoord, name=None)

Observation
    - target : FixedTarget
    - times : Time
    - gstimes :  Longitude (hourangle)
    - band : str
    - wavelength : u.Quantity
    - frequency : u.Quantity
    - datarate : u.Quantity
    - subbands : int
    - channels : int
    - polarizations : int
    - inttime : u.Quantity
    - ontarget_fraction : float
    - bandwidth : u.Quantity
    - bitsampling : u.Quantity
    - stations : Stations
    + __init__(target=None, times=None, band=None, datarate=None,subbands=None,
               channels=None, polarizations=None, inttime=None, ontarget=1.0,
               stations=None, bits=2)
    + elevations() --> dict[codename]: list
    + altaz() --> dict[codename]: list
    + is_visible() --> dict[codename]: list
    + bandwidth_smearing() --> u.Quantity
    + time_smearing() --> u.Quantity
    + datasize() --> u.Quantity
    + thermal_noise() --> u.Quantity
    + get_uv() --> dict[baseline-codename]: np.array (lambda units)

    <!-- + \_get_baseline_numbers() -->
    <!-- + \_get_baseline_number(ant1, ant2)  -> int -->



### plots.py


### functions.py


+ get_stations_from_file(filename) --> dict
+ stations_with_band(networks : dict/Stations, band) --> Stations


+ get_bandwidth_smearing()
+ get_time_smearing()



### freqsetups.py

Band()
    - name
    - frequency

Bands([fromfile])
    - names
    - frequencies




## For the doc:

- See the summary part, modal popover... Like accordion but only text from Bootstrap. Useful for different sections.
- There is a `address` style in Bootstrap, probably to implement in doc.


## To Do list

- Strike antenna name text if does not have the given frequency.
- Put some "loading" times.
- Add GST times.
- In sensitivity: highlight antennas that cannot observe the source.
- Add some pop up (`card` on Bootstrap) on each antenna with info: picture?, codename, country, etc.
- observation.Observation.thermal_noise() can be highly optimized.
- Help buttons and explanatory dialog.
- Arecibo limits. Currently not shown.
- Add a per-baseline basis sensitivity? (maybe as a roll-over?)
- Also, something about the largest angular scale you are sensitive to.






