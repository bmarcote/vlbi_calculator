# Observation

The `Observation` class is the main entry point for planning VLBI observations: it computes
source visibility, elevations, thermal noise (rms), baseline sensitivities, (u,v) coverage,
and synthesized beam size for a given network of stations and a set of scans.

An `Observation` is built from a `band`, a `Stations` network, and a `dict[str, ScanBlock]` of
scans to observe (all three are required at construction time — the object is not meant to be
built incrementally field-by-field, though `times`, `datarate`, and correlator settings can be
changed afterwards).

## Overview

```python
from vlbiplanobs.stations import Stations
from vlbiplanobs.sources import Source, Scan, ScanBlock
from vlbiplanobs.observation import Observation
from astropy.time import Time
from astropy import units as u

stations = Stations().filter_networks('EVN')
target = Source('M87')
scans = {'M87': ScanBlock([Scan(target)])}

obs = Observation(band='6cm', stations=stations, scans=scans,
                  times=Time('2025-03-15 08:00') + range(0, 481, 5)*u.min,
                  datarate=1024*u.Mbit/u.s)
```

## Key Properties

### Observing Setup

- `band` - Observing wavelength (e.g., '6cm', '18cm').
- `times` - Array of observation times (`astropy.time.Time`); setting it also fixes `duration` and `fixed_time`.
- `duration` - Total observation duration; used when `times` has not been fixed.
- `stations` - `Stations` network participating in the observation.
- `datarate` - Data rate per antenna (`astropy.units.Quantity`).
- `subbands`, `channels`, `polarizations`, `inttime`, `bitsampling`, `ontarget_fraction` - Correlator/observing setup parameters.
- `bandwidth` - Total recorded bandwidth, derived from `datarate`, `polarizations`, and `bitsampling`.

### Sources

- `scans` - Dictionary of `ScanBlock` objects keyed by block name.
- `sources(source_type=None)` - Returns all `Source` objects in the observation, optionally filtered by `SourceType`.
- `sourcenames` - Names of all sources in the observation.
- `sources_in_block(block_name, source_type=None)` / `sourcenames_in_block(...)` - Same, restricted to one scan block.

### Visibility

- `elevations()` - Elevation of each source for each station, at every observing time.
- `altaz()` - Full alt/az coordinates per source/station.
- `is_observable(times=None)` - Per scan-block, per-station boolean visibility (sources within a block are ANDed together).
- `per_source_observable()` - Same as above but per individual source (not ANDed within a block).
- `can_be_observed()` - Whether each scan block can be observed by each station at least some of the time.
- `is_observable_by_network(min_stations=2)` - Whether each scan block is observable by at least `min_stations` stations simultaneously.
- `is_always_observable()` - Whether each source is visible for the *entire* observation, per station.
- `when_is_observable(...)` - Suggests time windows when a source becomes observable.
- `sun_constraint(times=None)` / `sun_constraint_per_source(times=None)` - Minimum Sun separation checks.

### Results

- `thermal_noise()` - Expected rms thermal noise (baseline-sensitivity-weighted); returns a single `Quantity`, a per-source dict, or `None` if the band/data rate is not usable.
- `baseline_sensitivity(antenna1=None, antenna2=None)` - rms noise for a single baseline, or all baselines if antennas are not specified.
- `synthesized_beam()` - Estimated synthesized beam (`bmaj`, `bmin`, `pa`) per source, from an ellipse fit to the (u,v) coverage.
- `longest_baseline()` / `shortest_baseline()` - Longest/shortest projected baseline per source.
- `bandwidth_smearing()` / `time_smearing()` - Expected smearing effects.
- `datasize()` - Estimated total recorded data volume.
- `get_uv_data()` / `get_uv_values()` - Per-baseline and combined (u,v) points per source.
- `print_obs_times(date_format='%d %B %Y')` - Human-readable summary of the observation's time range.

### Scheduling

- `schedule_file(experiment_code='EXCODE', pi_name='PI Name', pi_email='pi@example.com', pi_institute='Institute', setup_file=None, comments='')` -
  Generates a pySCHED-compatible `.key` file as a string. This is the real entry point for
  turning an `Observation` into a schedule: internally it builds a
  `vlbiplanobs.scheduler.ObservationScheduler` for `self`, runs `.schedule()`, and uses the
  resulting scan placement to fill in the key-file template.

!!! warning "Not implemented"
    - `Observation.scheduler()` is a stub that unconditionally raises `NotImplementedError`.
      Do not use it. Use `schedule_file()` (or drive `vlbiplanobs.scheduler.ObservationScheduler`
      directly, see [Scheduler](scheduler.md)) instead.
    - `Observation.get_dirtymap()` also unconditionally raises `NotImplementedError` — the
      dirty-map/imaging code after the raise is dead and not usable. Do not document or rely on it.

## Examples

### Basic Observation

```python
from vlbiplanobs.stations import Stations
from vlbiplanobs.sources import Source, Scan, ScanBlock
from vlbiplanobs.observation import Observation
from astropy.time import Time
from astropy import units as u

stations = Stations().filter_networks('EVN')
scans = {'M87': ScanBlock([Scan(Source('M87'))])}

obs = Observation(band='6cm', stations=stations, scans=scans,
                  times=Time('2025-03-15 08:00') + range(0, 481, 5)*u.min,
                  datarate=1024*u.Mbit/u.s)

# Check visibility
visibility = obs.is_observable()
elevations = obs.elevations()
rms = obs.thermal_noise()
```

### Generate Schedule File

```python
key_content = obs.schedule_file(
    experiment_code='EG123A',
    pi_name='Your Name',
    pi_email='you@example.com',
)

with open('eg123a.key', 'w') as f:
    f.write(key_content)
```

### Multiple Targets

```python
scans = {
    'M87': ScanBlock([Scan(Source('M87'))]),
    '3C273': ScanBlock([Scan(Source('3C273'))]),
    'CygnusA': ScanBlock([Scan(Source('Cygnus A'))]),
}
obs.scans = scans

for name, block in obs.scans.items():
    print(f"{name}: {block.sourcenames()}")
```
