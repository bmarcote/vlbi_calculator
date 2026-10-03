# Stations

Classes for managing VLBI stations and networks.

## Station

Represents a single radio telescope.

### Key Properties

- `name` - Full station name
- `codename` - Two-letter station code (e.g., 'Ef', 'Wb')
- `networks` - Networks this station belongs to
- `location` - Geographic coordinates
- `bands` - Available observing bands
- `group` - Group identifier for multi-configuration antennas (e.g. `'VLA'`). `None` for standalone stations. Used by the GUI to collapse configurations into a single chip.
- `horizon` - Azimuth-dependent local horizon as a `(az_deg, el_deg)` tuple of NumPy arrays. `None` if no horizon data is available. Used to exclude elevations below terrain/structure blockage at a given azimuth.
- `horizon_min_elevation(az)` - Returns the minimum observable elevation (degrees) at a given azimuth. Returns `0` if no horizon is defined.

### Example

```python
from vlbiplanobs.stations import Stations

all_stations = Stations()          # reads the default station catalog
effelsberg = all_stations['Ef']

print(f"Name: {effelsberg.name}")
print(f"Location: {effelsberg.location}")
print(f"Bands: {effelsberg.bands}")
print(f"SEFD at 6 cm: {effelsberg.sefd('6cm')}")
```

`Stations(filename='my_stations.inp')` reads a custom catalog instead; a missing or invalid file raises an error.

## Stations

Collection of Station objects with filtering capabilities.

### Key Methods

- `filter_networks(networks, only_defaults=False)` - Stations that belong to the given network(s)
  (with `only_defaults=True`, only the default stations of each network).
- `stations_with_band(band)` - Iterate over the stations that can observe a band.
- `add_station(station)` / `remove_station(station)` - Modify the collection.
- `station_codenames`, `station_names`, `number_of_stations`, `observing_bands` - Properties.

### Example

```python
from vlbiplanobs.stations import Stations

# Load all stations
stations = Stations()

# Filter by network
evn = stations.filter_networks('EVN', only_defaults=True)
print(f"EVN stations: {evn.station_codenames}")

# Filter by band
cm6_capable = [s.codename for s in stations.stations_with_band('6cm')]

# Combine networks
combined = stations.filter_networks(['EVN', 'eMERLIN'])
```

The default stations of each network (as used by `planobs -n`) are available in `vlbiplanobs.NETWORKS`.

## Available Networks

| Network | Description |
|---------|-------------|
| EVN | European VLBI Network |
| eMERLIN | Enhanced Multi Element Remotely Linked Interferometer Network |
| VLBA | Very Long Baseline Array |
| HSA | High Sensitivity Array |
| LBA | Australian Long Baseline Array |
| KVN | Korean VLBI Network |
| VERA | VLBI Exploration of Radio Astrometry |
| KaVA | KVN and VERA Array |
| EAVN | East Asian VLBI Network |
| GMVA | Global mm-VLBI Array |
| EHT | Event Horizon Telescope |
| SKA-AA*, SKA-AA4 | SKA-Mid phased-up configurations |

Run `planobs --list-networks` for the current list and default stations.

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
