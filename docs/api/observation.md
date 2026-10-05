# Observation

The `Observation` class (`vlbiplanobs.observation.Observation`) holds the full definition of a VLBI observation
(band, stations, times, scan blocks, correlator setup) and computes all derived quantities.

## Creating an Observation

The recommended way is `vlbiplanobs.cli.main`, which resolves the targets and networks exactly as the CLI does:

```python
from astropy.time import Time
from astropy import units as u
from vlbiplanobs import cli

obs = cli.main(band='6cm', networks=['EVN'], stations=['Ar'], targets=['J1230+1223'],
               start_time=Time('2025-03-15 20:00', scale='utc'), duration=8*u.h)
```

Main keyword arguments of `cli.main`: `band`, `networks`, `stations`, `station_catalog`, `src_catalog`, `targets`,
`start_time` (UTC), `duration` (> 0 and ≤ 96 h), `datarate`, `ontarget`, `subbands`, `channels`, `polarizations`,
`inttime`, `phasecal_names`, `check_source_names`, `fringefinder_spec`, `polcal`. Without `start_time`/`duration`, the
observation covers a full GST day so that the best observing window can be searched.

## Key Properties

### Observing Setup

- `band` - Observing band (e.g. `'6cm'`, `'18cm'`); `wavelength` and `frequency` are derived from it.
- `times` - Observing times (`astropy.time.Time` array); `duration`, `gstimes` and `fixed_time` are derived from it.
- `stations` - `Stations` participating in the observation.
- `datarate` - Data rate (Mbit/s). Per-station rates are available with `station_datarate(codename)` and
  `station_datarates`; the shared `Station` objects are never modified per observation.
- `subbands`, `channels`, `polarizations`, `inttime`, `bandwidth`, `bitsampling`, `ontarget_fraction`.

### Sources

- `scans` - Dictionary of `ScanBlock` objects (one per target block).
- `sources(source_type=None)` - All sources in the observation (optionally filtered by `SourceType`).
- `sourcenames` - Names of all sources (a property).

### Results

- `thermal_noise()` - Expected image RMS noise per target block.
- `synthesized_beam()` - Approximate synthesized beam (`bmaj`, `bmin`, `pa`) per target block. Blocks without
  uv coverage are omitted, so use `.get(name)`.
- `elevations()` / `altaz()` - Source elevations / AltAz coordinates per station.
- `is_observable()` - Visibility of each source per station at each time.
- `when_is_observable(min_stations=3)` - Time ranges when at least `min_stations` stations can observe.
- `longest_baseline()`, `shortest_baseline()`, `bandwidth_smearing()`, `time_smearing()`, `datasize()`.
- `sun_constraint()` - Minimum Sun separation for each source (None if no constraint applies).
- `schedule_file(...)` - Text of a SCHED `.key` file (see [Scheduler](scheduler.md)).

## Examples

### Check visibility

```python
visibility = obs.is_observable()     # {source: {station: [bool, ...]}}
elevations = obs.elevations()        # {source: {station: Latitude array}}
windows = obs.when_is_observable()   # {source: [Time ranges]}
```

### Generate Schedule File

```python
key_content = obs.schedule_file(experiment_code='EG123A', pi_name='Your Name', email='you@example.com')

with open('eg123a.key', 'w') as f:
    f.write(key_content)
```

### Multiple Targets

```python
obs = cli.main(band='18cm', networks=['EVN'], targets=['J1230+1223', 'J1229+0203'],
               start_time=Time('2025-03-15 20:00', scale='utc'), duration=8*u.h)

for name, block in obs.scans.items():
    print(f"{name}: {block.sourcenames()}")
```
