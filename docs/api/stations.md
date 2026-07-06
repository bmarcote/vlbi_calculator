# Stations

Classes for managing VLBI stations and networks.

## Station

Represents a single radio telescope (antenna).

### Key Properties

- `name` - Full station name
- `codename` - Two-letter station code (e.g., 'Ef', 'Wb')
- `fullname` - Expanded name (e.g. 'Karl G. Jansky Very Large Array' for 'VLA')
- `networks` - Networks this station belongs to
- `country` - Country where the station is located
- `diameter` - Free-format string describing the antenna size (e.g. `'25 x 20 m'` for connected interferometers)
- `location` - Geographic coordinates (`astropy.coordinates.EarthLocation`)
- `bands` - Available observing bands
- `sefds` - Dict of system equivalent flux density (SEFD) per band
- `mount` - The `Mount` (mount type, axis limits, slewing speed/acceleration)
- `real_time` - Whether the station can join real-time (e.g. e-EVN) observations
- `decommissioned` - Whether the station is decommissioned
- `max_datarate` / `datarate` - Maximum and assumed recording data rates
- `sched_name` - Station name as used in the SCHED software catalog
- `group` - Group identifier for multi-configuration antennas (e.g. `'vla'`). `None` for standalone stations. Used by the GUI to collapse configurations into a single chip.
- `horizon` - Azimuth-dependent local horizon as a `(az_deg, el_deg)` tuple of NumPy arrays. `None` if no horizon data is available. Used to exclude elevations below terrain/structure blockage at a given azimuth.

### Key Methods

- `horizon_min_elevation(az)` - Minimum observable elevation (degrees) at a given azimuth. Returns `0` if no horizon is defined.
- `has_band(band)` - Whether the station can observe at a given band.
- `sefd(band)` - SEFD (in Jy) at a given band.
- `elevation(times, target)` / `altaz(times, target)` / `hour_angle(times, target)` - Compute elevation, Alt/Az, or hour angle of a source as seen from the station.
- `is_observable(times, target)` - Per-timestamp visibility given the station's constraints.
- `is_ever_observable(times, target)` / `is_always_observable(times, target)` - Whether a source is visible at any/every point in a time range.
- `slewing_time(coords1, coords2, time)` - Expected slewing time between two Alt/Az positions given the mount's speed and acceleration.

### Example

```python
from vlbiplanobs.stations import Stations

all_stations = Stations()
effelsberg = all_stations['Ef']

print(f"Name: {effelsberg.name}")
print(f"Location: {effelsberg.location}")
print(f"Bands: {effelsberg.bands}")
```

## Stations

Collection of `Station` objects representing a network of antennas, with filtering capabilities.

### Key Properties

- `name` - Full (expanded) name of the network
- `stations` - List of all `Station` objects in the network
- `number_of_stations` - Number of stations
- `station_names` / `station_codenames` - Names / codenames of all stations
- `observing_bands` - All bands observable by at least one station in the network

### Key Methods

- `filter_networks(networks, only_defaults=False)` - Get a new `Stations` containing only the antennas belonging to the given network(s) (e.g. `'EVN'`, or a list/comma-separated string of several).
- `filter_band(band)` - Get a new `Stations` containing only stations that can observe a given band.
- `filter_antennas(codenames)` - Get a new `Stations` containing only the stations with the given codenames.
- `add_station(station)` / `remove_station(station)` - Add or remove a `Station` (or iterable of `Station`) from the network.
- `stations_with_band(band)` - Generator yielding the stations that can observe at a given band.
- `has_band(band)` / `max_datarate(band)` - Check band availability / maximum data rate for the network.
- `Stations.get_networks_from_configfile()` (static) - Returns a `dict[str, Stations]` with all networks defined in the default (or given) network catalog file.
- `Stations.get_network_full_name(network)` (static) - Returns the expanded name of a network nickname.

Two `Stations` objects can also be combined with `+` (e.g. `evn + emerlin`), which merges their stations and intersects their observing bands.

### Example

```python
from vlbiplanobs.stations import Stations

# Load all stations
stations = Stations()

# Filter by network
evn = stations.filter_networks('EVN')
print(f"EVN stations: {evn.station_codenames}")

# Filter by band
cm6_capable = stations.filter_band('6cm')

# Combine networks
combined = stations.filter_networks(['EVN', 'eMERLIN'])
```

## Available Networks

| Network | Description |
|---------|-------------|
| EVN | European VLBI Network |
| eMERLIN | Enhanced Multi Element Remotely Linked Interferometer Network |
| VLBA | Very Long Baseline Array |
| LBA | Australian Long Baseline Array |
| KVN | Korean VLBI Network |
| GMVA | Global mm-VLBI Array |

## Station Codes

Common EVN station codes:

| Code | Station |
|------|---------|
| Ef | Effelsberg (Germany) |
| Wb | Westerbork (Netherlands) |
| Jb | Jodrell Bank (UK) |
| On | Onsala (Sweden) |
| Mc | Medicina (Italy) |
| Nt | Noto (Italy) |
| Tr | Torun (Poland) |
| Ys | Yebes (Spain) |
| Hh | Hartebeesthoek (South Africa) |
| Sh | Shanghai (China) |
