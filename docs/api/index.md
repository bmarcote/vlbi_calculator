# API Reference

Python API of the `vlbiplanobs` package. The same functionality as the `planobs` CLI is available from Python.

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

The simplest way to build a fully configured observation is `vlbiplanobs.cli.main` (also exported as
`vlbiplanobs.VLBIObs`), which takes the same inputs as the `planobs` CLI and returns an `Observation`:

```python
from astropy.time import Time
from astropy import units as u
from vlbiplanobs import cli

obs = cli.main(band='6cm', networks=['EVN'], targets=['J1230+1223'],
               start_time=Time('2025-03-15 20:00', scale='utc'), duration=8*u.h)

print(obs.thermal_noise())        # {'J1230+1223': <Quantity ... Jy / beam>}
print(obs.synthesized_beam())     # {'J1230+1223': {'bmaj': ..., 'bmin': ..., 'pa': ...}}

# Generate a SCHED .key file
key_content = obs.schedule_file(experiment_code='EG123A')
```

`start_time` must be a UTC `Time`, and `duration` must be > 0 and ≤ 96 h. A `ValueError` is raised otherwise,
or when no network/stations are given or a target cannot be resolved.

## Module Index

| Module | Description |
|--------|-------------|
| `vlbiplanobs.cli` | `planobs` command line; `main()` builds an `Observation` from CLI-like inputs |
| `vlbiplanobs.observation` | Main `Observation` class |
| `vlbiplanobs.stations` | `Station` and `Stations` classes |
| `vlbiplanobs.sources` | `Source`, `Scan`, `ScanBlock`, `SourceCatalog` classes |
| `vlbiplanobs.scheduler` | `ObservationScheduler` for pySCHED output |
| `vlbiplanobs.calibrators` | Radio Fundamental Catalog (`RFCCatalog`), fringe finder and phase calibrator searches |
| `vlbiplanobs.nme` | Network Monitoring Experiment planning (`plan_nme`, `generate_nme_key_file`) |
| `vlbiplanobs.freqsetups` | Frequency setup definitions |
| `vlbiplanobs.report` | PDF/TXT/Markdown/JSON observation reports |
