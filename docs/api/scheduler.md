# Scheduler

Observation scheduling for generating pySCHED-compatible files. This is the real scheduling
engine — `Observation.scheduler()` is only a stub (raises `NotImplementedError`); actual
scheduling happens through `ObservationScheduler`, either directly or via
`Observation.schedule_file()`, which builds one internally.

## ObservationScheduler

Main class for scheduling VLBI observations. Arranges fringe finders, polarisation calibrators,
eMERLIN 3C286, and science targets across the observation's time range, then generates SCHED
`.key` files (using `group N rep R` syntax and Jb1 source-change mitigation where relevant).

```python
ObservationScheduler(observation, min_antennas=2, require_all_antennas=False,
                     fringefinder_spec=None, polcal=False)
```

- `observation` - The `Observation` containing the scan blocks to schedule.
- `min_antennas` - Minimum antennas required for a valid time slot.
- `require_all_antennas` - If `True`, only schedule when all antennas can observe.
- `fringefinder_spec` - Source names, or `['N']` to auto-select *N* fringe-finder sources.
- `polcal` - Whether polarisation calibration scans are required.

If the `ortools` package is installed, scheduling prefers a CP-SAT constraint-programming
layout (`schedule()` falls back to the heuristic scheduler if CP-SAT finds no solution or
`ortools` is unavailable).

### Scheduling Rules

**Fringe Finders:**

- 2x 5-minute scans at observation start
- 1x 5-minute scan at observation end
- For observations > 3 hours: additional scans every ~2 hours

**Science Blocks:**

- For a single science block, scheduled in the gaps between fringe finders.
- For multiple science blocks, a dedicated multi-source mode equalises time per source,
  picks optimal visibility windows, minimises slewing, and interleaves fringe-finder scans
  between source blocks.

### Key Methods

- `schedule()` - Runs the scheduler and returns a `dict[str, ScanBlock]` mapping `'NNN_label'`
  ordered block names to the scheduled `ScanBlock`s. Also populates the internal scheduled-block
  list used by `get_scheduled_blocks()` and `generate_key_file()`.
- `get_scheduled_blocks()` - Returns the scheduled blocks (`list[ScheduledScanBlock]`) in
  chronological order. Must be called after `schedule()`.
- `generate_key_file(experiment_code='EXCODE', pi_name='PI Name', pi_email='pi@example.com', pi_institute='Institute', setup_file=None, comments='')` -
  Builds the full SCHED `.key` file content as a string from the current schedule.

### Example

```python
from vlbiplanobs.scheduler import ObservationScheduler

# Create scheduler
scheduler = ObservationScheduler(observation, min_antennas=3, require_all_antennas=False)

# Generate schedule
schedule = scheduler.schedule()

# Get detailed block info
for block in scheduler.get_scheduled_blocks():
    print(f"{block.name}: {block.start_time.iso} - {block.end_time.iso}")
    print(f"  Antennas: {block.n_antennas}, Elevation: {block.mean_elevation:.1f} deg")

# Write out a SCHED key file
key_content = scheduler.generate_key_file(experiment_code='EG123A', pi_name='Your Name')
with open('eg123a.key', 'w') as f:
    f.write(key_content)
```

## ScheduledScanBlock

Dataclass representing a scan block placed at a specific time, with timing and quality metadata.

### Properties

- `name` - Block identifier (e.g., `'FF_1'`, `'001_M87'`).
- `block` - Original `ScanBlock` being scheduled.
- `start_time` - Scheduled start time.
- `end_time` - Scheduled end time.
- `duration` - Computed property: `end_time - start_time`, in minutes.
- `scans` - Expanded list of individual `Scan`s filling this time slot.
- `n_antennas` - Number of antennas that can observe during this block.
- `mean_elevation` - Mean source elevation across observing antennas (degrees).

## Integration with Observation

The scheduler is automatically used when generating schedule files through `Observation`:

```python
from vlbiplanobs.observation import Observation

obs = Observation(band='6cm', stations=stations, scans=scans, times=times, datarate=datarate)

# This internally builds an ObservationScheduler(obs), calls .schedule(), and uses the result
key_content = obs.schedule_file(experiment_code='EG123A', pi_name='Your Name')
```

## Manual Scheduling

For more control over the scheduling process, use `ObservationScheduler` directly instead of
going through `Observation.schedule_file()`:

```python
from vlbiplanobs.scheduler import ObservationScheduler

scheduler = ObservationScheduler(obs, min_antennas=4, require_all_antennas=True)
schedule_dict = scheduler.schedule()  # dict[str, ScanBlock]

# Access scheduled blocks
blocks = scheduler.get_scheduled_blocks()

for block in blocks:
    print(f"{block.name}:")
    print(f"  Time: {block.start_time.iso} to {block.end_time.iso}")
    print(f"  Duration: {block.duration}")
    print(f"  Antennas: {block.n_antennas}")
    print(f"  Scans: {len(block.scans)}")

key_content = scheduler.generate_key_file(experiment_code='EG123A', pi_name='Your Name')
```
