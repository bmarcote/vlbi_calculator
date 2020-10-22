
#  EVN Source Elevation App


## To Do

[ ] Update the 'install_requires' in setup.py (line 54).
[ ] Move the app.py into bin/ folder, so it is directly called from terminal.
    https://python-packaging.readthedocs.io/en/latest/command-line-scripts.html





## Program structure

- app.py
- README.md
- docs/
- vlbiplanobs/
    - stations.py
    - sources.py
    - plots.py
    - functions.py
    - graphical_elements.py
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

+ get_doc_text()



### graphical_elements.py

+ tooltip(message, idname, trigger='?', placement='right', \*\*kwargs)
+ tooltip_card(a_card, idname, trigger='?', placement='right', \*\*kwargs)
+ create_accordion_card(title, text, id, is_open=True)
+ antenna_card(app, station)
+ antenna_cards(app, stations)
+ baseline_img(app, is_long=True)


### summary_cards.py




### stations.py

Station
    - observer : Observer
    - name : str
    - fullname : str
    - country : str
    - all_networks : str
    - codename : str
    - network : str
    - location : EarthLocation
    - bands : list
    - sefds : dict
    - min_elevation
    + __init__(name, codename, location, freqs_sefds, min_elevation=20, fullname=None,
all_networks=None, country=None)
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

SourceNotVisible(Exception)

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
    - duration : Time
    - band : str
    - wavelength : u.Quantity
    - frequency : u.Quantity
    - datarate : u.Quantity
    - subbands : int
    - channels : int
    - polarizations : int
    - inttime : u.Quantity
    - ontarget_fraction : float
    - ontarget_time : Time
    - bandwidth : u.Quantity
    - bitsampling : u.Quantity
    - stations : Stations
    + __init__(target=None, times=None, band=None, datarate=None,subbands=None,
               channels=None, polarizations=None, inttime=None, ontarget=1.0,
               stations=None, bits=2)
    + elevations() --> dict[codename]: list
    + altaz() --> dict[codename]: list
    + is_visible() --> dict[codename]: list
    + longest_baseline() --> (str, u.Quantity)
    + shortest_baseline() --> (str, u.Quantity)
    + bandwidth_smearing() --> u.Quantity
    + time_smearing() --> u.Quantity
    + datasize() --> u.Quantity
    + thermal_noise() --> u.Quantity
    + get_uv_baseline() --> dict[baseline-codename]: np.array (lambda units)
    + get_uv_array() --> np.array (lambda units)

    <!-- + \_get_baseline_numbers() -->
    <!-- + \_get_baseline_number(ant1, ant2)  -> int -->



### plots.py


### functions.py


+ get_stations_from_file(filename) --> dict
+ stations_with_band(networks : dict/Stations, band) --> Stations
+ print_obs_times(obs, dateformat='%d &m &Y')




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
X Add GST times.
X In sensitivity: highlight antennas that cannot observe the source.
- Add some pop up (`card` on Bootstrap) on each antenna with info: picture?, codename, country, etc.
- observation.Observation.thermal_noise() can be highly optimized.
- Help buttons and explanatory dialog.
- Arecibo limits. Currently not shown.
- Add a per-baseline basis sensitivity? (maybe as a roll-over?)
- Also, something about the largest angular scale you are sensitive to.
- Add favicon.ico to the assets folder.






