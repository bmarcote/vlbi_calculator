# Quick Start

Get started with PlanObs in under 5 minutes.

## Your First Observation

### Using the CLI

Check when a source is visible with the EVN:

```bash
planobs -b 6cm -t 'M87' --network EVN
```

This will show:

- Source visibility for each antenna
- Expected sensitivity (RMS noise)
- Angular resolution

### Specify a Time

Plan an observation at a specific date:

```bash
planobs -b 6cm -t 'Cygnus A' --network EVN --starttime '2025-03-15 08:00' --duration 8
```

### Add More Stations

Combine networks and individual antennas:

```bash
planobs -b 6cm -t 'M87' --network EVN eMERLIN --stations Ar Ys
```

## Using the GUI

Launch the web interface:

```bash
planobs-server
```

Then open your browser to `http://localhost:8050`.

## Using Python

```python
from vlbiplanobs import Observation
from astropy.time import Time
from astropy import units as u

# Create observation
obs = Observation()
obs.band = '6cm'
obs.times = Time('2025-03-15 08:00') + (range(0, 481, 5) * u.min)

# Add EVN stations
obs.stations = obs.stations.filter_by_network('EVN')

# Add a target
obs.add_target('M87')

# View results
print(obs.sensitivity())
print(obs.angular_resolution())
```

## Generate a Schedule File

Create a pySCHED-compatible `.key` file:

```bash
planobs -b 6cm -t 'M87' --network EVN --starttime '2025-03-15 08:00' --duration 8 --sched my_observation
```

This creates `my_observation.key` ready for pySCHED.

## Next Steps

- [Tutorial](tutorial.md) – Complete usage guide
- [CLI Reference](cli.md) – All command-line options
- [API Reference](api/index.md) – Python API documentation
