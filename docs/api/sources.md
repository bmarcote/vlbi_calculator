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
```

### Example

```python
from vlbiplanobs.sources import Source, SourceType
from astropy.coordinates import SkyCoord

# Create from coordinates
src = Source(
    name='M87',
    coord=SkyCoord('12h30m49.4s', '+12d23m28s'),
    type=SourceType.TARGET
)

# Resolve by name (uses SIMBAD/NED)
src = Source.from_name('Cygnus A', type=SourceType.TARGET)
```

## Scan

A single observation of a source for a specific duration.

### Example

```python
from vlbiplanobs.sources import Scan
from astropy import units as u

scan = Scan(source=my_source, duration=10*u.min)
```

## ScanBlock

A collection of scans to be observed together.

### Key Methods

- `fill(duration)` - Expand scans to fill available time
- `sources(type)` - Get sources of a specific type
- `has(type)` - Check if block contains a source type

### Example

```python
from vlbiplanobs.sources import ScanBlock, Scan, Source, SourceType
from astropy import units as u

# Create scans
target = Scan(Source.from_name('M87', SourceType.TARGET), duration=10*u.min)
phasecal = Scan(Source.from_name('J1230+1223', SourceType.PHASECAL), duration=2*u.min)

# Create block
block = ScanBlock([phasecal, target])

# Check contents
print(f"Has target: {block.has(SourceType.TARGET)}")
print(f"Has fringe finder: {block.has(SourceType.FRINGEFINDER)}")

# Fill time slot
expanded_scans = block.fill(60*u.min)
```

## SourceType Enum

