# Sources

Classes for astronomical sources and observation scans.

## Source

Represents an astronomical source to observe.

### Source Types

```python
from vlbiplanobs.sources import SourceType

SourceType.TARGET        # Science target
SourceType.PHASECAL      # Phase calibrator
SourceType.FRINGEFINDER  # Fringe finder
SourceType.AMPLITUDECAL  # Amplitude calibrator
SourceType.CHECKSOURCE   # Check source
SourceType.POLCAL        # Polarization calibrator
SourceType.PULSAR        # Pulsar
SourceType.UNKNOWN       # Default when not given
```

### Example

```python
from vlbiplanobs.sources import Source, SourceType

# Create from coordinates
src = Source('M87', '12h30m49.4s +12d23m28s', source_type=SourceType.TARGET)

# Resolve by name (RFC catalog first, then SIMBAD/NED/VizieR online)
src = Source.source_from_name('Cygnus A', source_type=SourceType.TARGET)

# Parse a CLI-like string: name, coordinates, or 'name/coordinates'
cal = Source.source_from_str('MyCal/12h31m00s +12d00m00s', source_type=SourceType.PHASECAL)
```

RFC lookups are exact (case-insensitive J2000 or IVS name). Names sent to online resolvers are validated first
(`vlbiplanobs.sources.validate_source_name`) and the results are cached.

## Scan

A single observation of a source for a specific duration.

### Example

```python
from vlbiplanobs.sources import Scan
from astropy import units as u

scan = Scan(source=src, duration=10*u.min)
```

## ScanBlock

A collection of scans to be observed together.

### Key Methods

- `fill(max_duration)` - Expand the scans (repeating cycles) to fill the available time.
- `sources(source_type=None)` / `sourcenames(source_type=None)` - Sources in the block.
- `has(source_type)` - Check if the block contains a source type.
- `fractional_time()` - Fraction of time spent on each source.

Scans have a `duration` (default 10 min) and `every` (`-1`: every cycle; `N > 0`: every N cycles, e.g. check sources).

### Example

```python
from vlbiplanobs.sources import ScanBlock, Scan, Source, SourceType
from astropy import units as u

# Create scans
target = Scan(Source('M87', '12h30m49.4s +12d23m28s', source_type=SourceType.TARGET), duration=3.5*u.min)
phasecal = Scan(Source.source_from_name('J1230+1223', SourceType.PHASECAL), duration=1.5*u.min)

# Create block
block = ScanBlock([phasecal, target])

# Check contents
print(f"Has target: {block.has(SourceType.TARGET)}")
print(f"Has fringe finder: {block.has(SourceType.FRINGEFINDER)}")

# Fill time slot
expanded_scans = block.fill(60*u.min)
```

## SourceCatalog

`SourceCatalog(personal_catalog='my_sources.toml')` reads a personal TOML source catalog (the `-sc` CLI option),
with `blocknames`, `targets`, `pulsars`, `fringefinders`, `polcals`, `ampcals` and `sources()`.
