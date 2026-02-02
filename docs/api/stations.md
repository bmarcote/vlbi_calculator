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

Collection of Station objects with filtering capabilities.

### Key Methods

- `filter_by_network(network)` - Get stations in a network
- `filter_by_band(band)` - Get stations that can observe a band
- `filter_antennas(codenames)` - Get specific stations by code

### Example

```python
from vlbiplanobs.stations import Stations

# Load all stations
stations = Stations()

# Filter by network
evn = stations.filter_by_network('EVN')
print(f"EVN stations: {evn.station_codenames}")

# Filter by band
cm6_capable = stations.filter_by_band('6cm')

# Combine networks
combined = stations.filter_by_network(['EVN', 'eMERLIN'])
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
