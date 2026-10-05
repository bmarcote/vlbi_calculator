# Scheduler

Observation scheduling for generating pySCHED-compatible files (`vlbiplanobs.scheduler`).

## ObservationScheduler

Main class for scheduling VLBI observations. It works on its own copy of the observation scans, so the
`Observation` passed to it is never modified.

```python
ObservationScheduler(observation, min_antennas=2, require_all_antennas=False,
                     fringefinder_spec=None, polcal=False)
```

- `fringefinder_spec` - List of fringe finder names (or `name/coordinates`), or a single number (as a string,
  default `['2']`): how many fringe-finder scans to place on an auto-selected source. Same format as the
  `--fringefinders` CLI option.
- `polcal` - Add polarization calibrator scans (3C84, OQ208, DA193).

### Scheduling Rules

**Fringe Finders:**

- 5-minute scans, spread across the observation
- With a number N in `fringefinder_spec`: N scans on one auto-selected source (default 2)
- With source names: 1 scan up to 1.5 h, 2 up to 3 h, and one more per additional ~2 h; sources are used round-robin

**Added calibrators** (fringe finders, polcals, eMERLIN 3C286) are placed using their real visibility from the
stations. Polcals that are not visible by enough antennas are skipped with a warning.

**Science Blocks:**

- Scheduled in gaps between fringe finders
- Optimized for maximum antenna participation
- Secondary optimization for highest elevation

### Example

```python
from vlbiplanobs.scheduler import ObservationScheduler

# Create scheduler
scheduler = ObservationScheduler(obs, min_antennas=3, require_all_antennas=False)

# Generate schedule ({block name: ScanBlock})
schedule = scheduler.schedule()

# Get detailed block info
for block in scheduler.get_scheduled_blocks():
    print(f"{block.name}: {block.start_time.iso} - {block.end_time.iso}")
    print(f"  Antennas: {block.n_antennas}, Elevation: {block.mean_elevation:.1f}°")

# Write the SCHED .key file
key_content = scheduler.generate_key_file(experiment_code='EG123A', pi_name='Your Name',
                                          pi_email='you@example.com', setup_file='evn6cm-2Gbps-32MHz.set')
```

`generate_key_file` strips quotes and newlines from all user-provided strings (experiment code, PI fields,
comments, source names), so they cannot break the SCHED syntax. `template_path` selects a custom template.

## ScheduledScanBlock

Dataclass representing a scheduled scan block with timing metadata.

### Properties

- `name` - Block identifier (e.g., 'FF_1', 'Target_A')
- `block` - Original ScanBlock
- `start_time` - Scheduled start time
- `end_time` - Scheduled end time
- `duration` - Block duration
- `scans` - Expanded list of individual scans
- `n_antennas` - Number of participating antennas
- `mean_elevation` - Mean source elevation (degrees)

## Integration with Observation

`Observation.schedule_file()` produces a simple `.key` file directly from the observation:

```python
key_content = obs.schedule_file(experiment_code='EG123A', pi_name='Your Name')
```

The `planobs --sched` CLI option uses `ObservationScheduler.generate_key_file()`.

## Manual Scheduling

For more control over the scheduling process:

```python
scheduler = ObservationScheduler(obs, min_antennas=4, require_all_antennas=True)
schedule_dict = scheduler.schedule()

for block in scheduler.get_scheduled_blocks():
    print(f"{block.name}:")
    print(f"  Time: {block.start_time.iso} to {block.end_time.iso}")
    print(f"  Duration: {block.duration}")
    print(f"  Antennas: {block.n_antennas}")
    print(f"  Scans: {len(block.scans)}")
```
