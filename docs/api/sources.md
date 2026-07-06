# Sources

Classes for astronomical sources and observation scans.

## Source

Represents an astronomical source to observe. Subclasses `astroplan.FixedTarget`.

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
SourceType.UNKNOWN       # Default, unspecified type
```

### Key Properties

- `type` - The source's `SourceType`
- `flux` - Optional `SourceFlux` with estimated flux densities per band
- `notes` - Optional free-text notes
- `other_names` - List of alternative names for the source

### Key Methods

- `Source.source_from_name(name, source_type=SourceType.TARGET)` (classmethod) - Build a `Source` by resolving its coordinates from a name (RFC catalog first, then ICRS/online resolvers).
- `Source.source_from_str(spec, source_type=SourceType.TARGET)` (classmethod) - Build a `Source` from a name, a coordinate string, or a combined `'name/coordinates'` spec.
- `sun_separation(times)` - Angular separation to the Sun at given times.
- `sun_constraint(min_separation, times=None)` - Times when the Sun is too close to the source to observe.

### Example

```python
from vlbiplanobs.sources import Source, SourceType
from astropy.coordinates import SkyCoord

# Create from explicit coordinates
src = Source(
    name='M87',
    coordinates=SkyCoord('12h30m49.4s', '+12d23m28s'),
    source_type=SourceType.TARGET
)

# Resolve by name (RFC catalog, then SIMBAD/NED/VizieR)
src = Source.source_from_name('Cygnus A', source_type=SourceType.TARGET)
```

## SourceFlux / FluxMeasurement

`SourceFlux` stores per-band flux information for a `Source`, as a dict of `FluxMeasurement` (with
`resolved` = peak brightness on the longest baselines, and `unresolved` = total flux density on the
shortest baselines).

Key methods: `bands()`, `flux_density(band)`, `peak_flux(band)`, `has_band(band)`, `add_band(band, flux)`.

## Scan

A single pointing to a source for a given duration, optionally repeated only every `every` cycles
within a `ScanBlock`.

### Example

```python
from vlbiplanobs.sources import Scan
from astropy import units as u

scan = Scan(source=my_source, duration=10*u.min)
```

## ScanBlock

A collection of `Scan` objects to be observed together (e.g. target + phase calibrator + check source).

### Key Methods

- `fill(max_duration)` - Expand the scans to fill the available time, bracketing targets with phase
  calibrators and honoring any `every=N` scheduling.
- `sources(source_type=None)` / `sourcenames(source_type=None)` - Sources (or their names) of a specific type, or all if omitted.
- `has(source_type)` - Check if the block contains a given source type.
- `scan_with_sourcename(name)` - The `Scan` pointing to the source with the given name.
- `scans_with_sources(source_type)` - All scans for a given source type.
- `fractional_time()` - Estimated fraction of the block's time spent on each source.

### Example

```python
from vlbiplanobs.sources import ScanBlock, Scan, Source, SourceType
from astropy import units as u

# Create scans
target = Scan(Source.source_from_name('M87', SourceType.TARGET), duration=10*u.min)
phasecal = Scan(Source.source_from_name('J1230+1223', SourceType.PHASECAL), duration=2*u.min)

# Create block
block = ScanBlock([phasecal, target])

# Check contents
print(f"Has target: {block.has(SourceType.TARGET)}")
print(f"Has fringe finder: {block.has(SourceType.FRINGEFINDER)}")

# Fill time slot
expanded_scans = block.fill(60*u.min)
```

## SourceCatalog

Holds the `ScanBlock` objects (targets, pulsars) read from a personal TOML catalog file.

- `SourceCatalog(personal_catalog=None)` - Optionally reads a TOML catalog file immediately.
- `read_personal_catalog(path)` - Reads (or re-reads) a TOML catalog file, populating `targets` and,
  if present, `pulsars`.
- `targets` - Dict `{name: ScanBlock}` of target sources.
- `pulsars` - Dict `{name: ScanBlock}` of pulsar sources, or `None` if there are none.
- `blocks` / `blocknames` - All `ScanBlock` entries across block types, or just their names.
- `source_names(include_calibrators=False)` / `sources(include_calibrators=False)` - Names, or
  `{name: Source}`, of all sources in the catalog (targets only unless `include_calibrators=True`).

!!! note
    `SourceCatalog.ampcals`, `.fringefinders`, and `.polcals` are defined on the class but are
    currently unused: `read_personal_catalog` does not populate those catalog blocks, so these
    properties always return `None`. `read_rfc_catalog()` is also not implemented yet and raises
    `NotImplementedError` if called.

### Example

```python
from vlbiplanobs.sources import SourceCatalog

catalog = SourceCatalog('my_sources.toml')
print(catalog.source_names())
target_block = catalog.targets['M87']
```

## SourceNotVisible

Exception raised when a target source cannot be observed from any antenna in the network.
