# Release Notes

## Version 4.7.0 (Current)

### New Features

- **Grouped antenna chips in GUI** – Antennas with multiple configurations (e.g. VLA, MeerKAT, SKAO) are now collapsed into a single split-button chip. The chip label shows the group name; a **▼** dropdown lets you switch between configurations. Picking a configuration automatically selects the antenna for the observation. The active configuration is highlighted in blue in the dropdown.
- **Azimuth-dependent local horizon** – Station-specific terrain and structure blockage is now modelled for 25 stations (Effelsberg, VLBA, Urumqi/Nanshan, Robledo, and others) using pySCHED's horizon data. Elevations are correctly excluded when below the local horizon at a given azimuth.
- **Terminal elevation plot** – The CLI per-source visibility output now uses a `plotext` scatter plot with colour-coded elevation bands (red `< 10°`, yellow `10–20°`, green `20–40°`, cyan `40–60°`, blue `> 60°`) instead of hand-drawn character rows.
- **Phase calibrator and check-source selection** – `planobs observe` now accepts `--phasecal` and `--check-source` options for automatic or named calibrator selection.
- **Fringe finder improvements** – `planobs fringefinders` accepts `-n/--network` to filter by network; fringe finder scheduling uses round-robin distribution across multiple sources.
- **Polarisation calibration** – `--polcal` schedules 3C84/OQ208/DA193 blocks spread at 10 %/50 %/90 % of available time.
- **Name/coordinate source parsing** – All source arguments accept a `name/RA Dec` syntax (e.g. `'MySrc/12:30:49 +12:23:28'`) where user-supplied coordinates override any catalog lookup.

### Improvements

- Renamed Urumqi station display name to **Nanshan** (its current operational name).
- Added `group` field to station catalog entries for VLA (Y1, Y27), MeerKAT (Me1, Me), and SKAO (Sk1, Sk2, Sk4).
- Observability summary wording clarified: *"All antennas can observe the block simultaneously"* and *"Optimal visibility window"*.
- GUI chip alignment: grouped chips now match `dmc.Chip` height, font size, and weight exactly.
- `ortools >= 9.0` is now a hard dependency (was commented out).

### Bug Fixes

- `enable_antennas_with_band` callback crashed on initial page load when `band_index` was `None`.
- Network switches could add both grouped codenames (e.g. Y1 and Y27) simultaneously to the selection, violating the single-active-config invariant.
- Grouped chip dropdown highlight was reset on observation re-render; fixed by using CSS class patching instead of Mantine prop updates.

### API Changes

- `Station` gained an optional `group: Optional[str]` parameter and property for GUI grouping.
- `Station` gained an optional `horizon: Optional[tuple[np.ndarray, np.ndarray]]` parameter and `horizon_min_elevation(az)` method.
- `Source.parse_source_spec(spec)` static method added for `name/coord` parsing.
- `calibrators.select_phase_calibrator()` and `calibrators.select_check_source()` added.

---

## Version 4.6.4

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
