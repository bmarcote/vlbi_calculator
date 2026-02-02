# Release Notes

## Version 4.6.4 (Current)

### New Features

- **Observation Scheduler** – New `ObservationScheduler` class for automated scan scheduling
- **Schedule File Generation** – Generate pySCHED-compatible `.key` files via `--sched` CLI option
- **Fringe Finder Scheduling** – Automatic placement of fringe finders (2 at start, 1 at end, ~2h intervals)

### Improvements

- Optimized scheduler code with reduced memory footprint
- NumPy-style docstrings across all modules
- Complete MkDocs documentation

### API Changes

- Added `Observation.schedule_file()` method
- Added `ObservationScheduler` class in `vlbiplanobs.scheduler`
- Added `ScheduledScanBlock` dataclass

## Version 4.6.x

### Features

- Modern Dash-based GUI
- Real-time visibility calculations
- Interactive elevation plots
- Multiple network support

## Version 4.5.x

### Features

- Initial Python 3.11+ support
- CLI interface (`planobs` command)
- Source resolution via SIMBAD/NED/VizieR

## Migration Guide

### From 4.5.x to 4.6.x

No breaking changes. New scheduling features are additive.

### From EVN Calculator

PlanObs is a complete rewrite. Key differences:

1. **Local Installation** – Can run offline
2. **Python API** – Full programmatic access
3. **Multiple Outputs** – GUI, CLI, and schedule files
4. **Extended Networks** – Support for non-EVN arrays
