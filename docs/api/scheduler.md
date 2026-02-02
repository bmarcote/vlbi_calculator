# Scheduler

Observation scheduling for generating pySCHED-compatible files.

## ObservationScheduler

Main class for scheduling VLBI observations.

### Scheduling Rules

**Fringe Finders:**

- 2× 5-minute scans at observation start
- 1× 5-minute scan at observation end
- For observations > 3 hours: additional scans every ~2 hours

**Science Blocks:**

- Scheduled in gaps between fringe finders
- Optimized for maximum antenna participation
- Secondary optimization for highest elevation

### Example

```python
from vlbiplanobs.scheduler import ObservationScheduler

# Create scheduler
scheduler = ObservationScheduler(
    observation,
    min_antennas=3,
    require_all_antennas=False
)

# Generate schedule
schedule = scheduler.schedule()

# View results
scheduler.print_schedule()

# Get detailed block info
for block in scheduler.get_scheduled_blocks():
    print(f"{block.name}: {block.start_time.iso} - {block.end_time.iso}")
    print(f"  Antennas: {block.n_antennas}, Elevation: {block.mean_elevation:.1f}°")
```

## ScheduledScanBlock

Dataclass representing a scheduled scan block with timing metadata.

### Properties

- `name` - Block identifier (e.g., 'FF_start_1', 'Target_A')
- `block` - Original ScanBlock
- `start_time` - Scheduled start time
- `end_time` - Scheduled end time
- `duration` - Block duration
- `scans` - Expanded list of individual scans
- `n_antennas` - Number of participating antennas
- `mean_elevation` - Mean source elevation (degrees)

## Integration with Observation

The scheduler is automatically used when generating schedule files:

```python
from vlbiplanobs import Observation

obs = Observation()
# ... configure observation ...

# This internally uses ObservationScheduler
key_content = obs.schedule_file(
    experiment_code='EG123A',
    pi_name='Your Name'
)
```

## Manual Scheduling

For more control over the scheduling process:

```python
from vlbiplanobs.scheduler import ObservationScheduler

scheduler = ObservationScheduler(obs, min_antennas=4, require_all_antennas=True)
schedule_dict = scheduler.schedule()

# Access scheduled blocks
blocks = scheduler.get_scheduled_blocks()

for block in blocks:
    print(f"{block.name}:")
    print(f"  Time: {block.start_time.iso} to {block.end_time.iso}")
    print(f"  Duration: {block.duration}")
    print(f"  Antennas: {block.n_antennas}")
    print(f"  Scans: {len(block.scans)}")
```
