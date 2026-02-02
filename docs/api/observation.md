# Observation

The `Observation` class is the main entry point for planning VLBI observations.

## Overview

```python
from vlbiplanobs import Observation

obs = Observation()
obs.band = '6cm'
obs.add_target('M87')
```

## Class Reference

## Key Properties

### Observing Setup

- `band` - Observing wavelength (e.g., '6cm', '18cm')
- `times` - Array of observation times
- `stations` - Stations participating in the observation
- `datarate` - Data rate in Mbit/s

### Sources

- `scans` - Dictionary of scan blocks
- `sources()` - Get all sources in the observation

### Results

- `sensitivity()` - Calculate expected RMS noise
- `angular_resolution()` - Calculate synthesized beam size
- `is_observable()` - Check source visibility per station

## Examples

### Basic Observation

```python
from vlbiplanobs import Observation
from astropy.time import Time
from astropy import units as u

obs = Observation()
obs.band = '6cm'
obs.times = Time('2025-03-15 08:00') + (range(0, 481, 5) * u.min)
obs.stations = obs.stations.filter_by_network('EVN')
obs.add_target('M87')

# Check visibility
visibility = obs.is_observable()
elevations = obs.elevations()
```

### Generate Schedule File

```python
key_content = obs.schedule_file(
    experiment_code='EG123A',
    pi_name='Your Name',
    pi_email='you@example.com'
)

with open('eg123a.key', 'w') as f:
    f.write(key_content)
```

### Multiple Targets

```python
obs.add_target('M87')
obs.add_target('3C273')
obs.add_target('Cygnus A')

for name, block in obs.scans.items():
    print(f"{name}: {block.sources()}")
```
