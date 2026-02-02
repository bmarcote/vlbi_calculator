# API Reference

Complete Python API documentation for vlbiplanobs.

## Core Modules

### Observation

The main class for planning VLBI observations.

[Full Observation Documentation →](observation.md)

### Stations

Classes for managing VLBI stations and networks.

[Full Stations Documentation →](stations.md)

### Sources

Classes for astronomical sources and scan blocks.

[Full Sources Documentation →](sources.md)

### Scheduler

Observation scheduling for pySCHED output.

[Full Scheduler Documentation →](scheduler.md)

## Quick Example

```python
from vlbiplanobs import Observation
from astropy.time import Time
from astropy import units as u

# Create and configure observation
obs = Observation()
obs.band = '6cm'
obs.times = Time('2025-03-15 08:00') + (range(0, 481, 5) * u.min)
obs.stations = obs.stations.filter_by_network('EVN')
obs.add_target('M87')

# Get results
print(f"Sensitivity: {obs.sensitivity()}")
print(f"Resolution: {obs.angular_resolution()}")

# Generate schedule file
key_content = obs.schedule_file(experiment_code='EG123A')
```

## Module Index

| Module | Description |
|--------|-------------|
| `vlbiplanobs.observation` | Main Observation class |
| `vlbiplanobs.stations` | Station and Stations classes |
| `vlbiplanobs.sources` | Source, Scan, ScanBlock classes |
| `vlbiplanobs.scheduler` | ObservationScheduler for pySCHED output |
| `vlbiplanobs.freqsetups` | Frequency setup definitions |
| `vlbiplanobs.rfc` | Radio Fundamental Catalog interface |
