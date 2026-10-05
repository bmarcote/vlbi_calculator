# Release Notes

## Version 5.1.0 (Current)

*Released 2026-10-03.* Security, correctness, and performance release. It also unifies the CLI option names across
all subcommands (old spellings keep working with a warning, see the migration table below).

### New Features

- **Unified CLI option names** – The same concept now uses the same flag in every subcommand: `-e/--epoch`
  (start time), `-t/--target`, `-n/--network`, `-s/--stations`, `-b/--band`, `-d/--duration`, `-l/--max-lines`,
  `-sc/--source-catalog`, `--station-catalog`, `--json`. See the [CLI overview](cli.md#option-names).
- `planobs phasecals` accepts `-sc/--source-catalog`: the target is first looked up in your personal source catalog.
- `planobs observe -d/--duration` accepts decimal hours (e.g. `-d 1.5`).
- **Polaris export allowlist** – The new `PLANOBS_POLARIS_ORIGINS` environment variable restricts which web origins
  can receive the setup via **Export to Polaris** (see [Web Server](mode-server.md#environment-variables)).

### Breaking or Deprecated CLI Changes

| Subcommand | Old | New | Status of the old spelling |
|------------|-----|-----|----------------------------|
| observe | `-t1`, `--starttime` | `-e`, `--epoch` | Deprecated (works, prints a warning) |
| observe | `--targets` | `-t`, `--target` | Deprecated (works, prints a warning) |
| fringefinders | `--starttime` | `-e`, `--epoch` | Deprecated (works, prints a warning) |
| fringefinders | `-t` | `-e`, `--epoch` | **Removed** (error: `-t` means `--target` elsewhere) |
| phasecals | `--n-sources` | `-l`, `--max-lines` | Deprecated (works, prints a warning) |
| phasecals | `-n` | `-l`, `--max-lines` | **Removed** (error: `-n` means `--network` elsewhere) |
| phasecals | `--catalog-file` | `--rfc-catalog` | Deprecated (works, prints a warning) |

Deprecated spellings are hidden from `-h` and will be removed in a future version. Update your scripts, e.g.:

```bash
# Before
planobs -b 6cm --targets M87 -n EVN -t1 '2025-03-15 08:00' -d 8
planobs fringefinders -s Ef Hh Mc Tr -t '2025-03-15 08:00' -d 8
planobs phasecals M87 -n 10 --catalog-file my_rfc.txt
# Now
planobs -b 6cm -t M87 -n EVN -e '2025-03-15 08:00' -d 8
planobs fringefinders -s Ef Hh Mc Tr -e '2025-03-15 08:00' -d 8
planobs phasecals M87 -l 10 --rfc-catalog my_rfc.txt
```

The `planobs fringefinders` and `planobs phasecals` entry points in `vlbiplanobs.calibrators` now reuse the same
parsers as `planobs`, so they accept exactly the same options.

### Security

- **Server-side validation of all GUI inputs** – Values from the GUI, shared `?config=` links, the browser's local
  storage, and uploaded files are validated on the server: duration ≤ 50 h, ≤ 20 targets, names ≤ 80 characters,
  uploads ≤ 50 kB, and only known setup values (band, data rate, subbands, channels, polarizations, integration time).
- **No cross-user data leaks** – Per-station data rates are stored per observation; the shared station catalog is no
  longer modified, so one user's setup cannot affect another user's results.
- **SCHED key file hardening** – Quotes and newlines are stripped from user-provided strings, placeholders are filled
  in a single pass, and `--sched` experiment codes may only contain letters, digits, and underscores.
- Source names are validated before online name resolution, and the results are cached.
- CLI output escapes user and catalog strings, so they cannot inject Rich markup.
- **Export to Polaris** only posts the setup to the exact opener origin (optionally restricted by
  `PLANOBS_POLARIS_ORIGINS`), never to `*`.
- **Mixpanel analytics removed** – PlanObs no longer uses any tracker.

### Bug Fixes

- RFC catalog lookups are exact (case-insensitive J2000 or IVS name); a substring match could return the wrong source.
- The scheduler no longer modifies the `Observation` it schedules, and added calibrators (eMERLIN 3C286, polarization
  calibrators, fringe finders) are placed using their real visibility. Generated `.key` files may therefore differ
  from previous versions; polarization calibrators that cannot be placed are skipped with a warning.
- Automatic phase calibrator / check source selection never picks the target itself under another name (e.g.
  RFC `J1230+1223` for `M87`): candidates within 5″ of the target or phase calibrator are skipped.
- The eMERLIN 3C286 flux-scale scan is no longer added to EVN-only observations. Only the stations Cm, Da, De, Kn and
  Pi trigger it; Jodrell Bank (Jb1, Jb2) observes regularly within the EVN and does not.
- The CLI time grid no longer duplicates the last time sample.
- The spurious `RuntimeWarning: invalid value encountered in do_format` printed when showing coordinates (an
  astropy/numpy ≥ 2.4 incompatibility) is silenced.
- Custom network and station catalogs (`--station-catalog`) are actually read, or fail with a clear error.
- Fixed eMERLIN-only arrays (per-band maximum data rates), sources that are never visible, visibility windows across
  midnight, scan blocks with check sources but without phase calibrator (`every=N`), TOML scans without a duration,
  and catalogs that only contain pulsars.
- AstroGeo links are correct for sources with −1° < Dec < 0°.
- GUI antenna highlighting uses exact matching (e.g. `Me` no longer highlights `Me1`).
- The GUI maximum-duration message now says 50 h, matching the actual limit.

### Performance

- Visibility and elevation computations are vectorized with ERFA (~17× faster).
- The RFC catalog is parsed once per process (2.2 s → 0.13 s).
- Scheduling is ~3–5× faster.
- The GUI sends the uv-coverage plot once per tab (~50 % smaller payload) and highlights antennas in the browser,
  without a server round-trip.

### Packaging

- Python 3.12+ is required.
- Added `pyerfa`; dropped the unused `six`, `matplotlib`, `types-PyYAML`, and `Cython` dependencies.
- `kaleido` is kept below 1 (kaleido ≥ 1 needs a Chrome binary to export figures to PDF), excluding 0.2.1.post1
  (no x86-64 wheels; `uv` could pick it on a fresh install).
- `plotext` is kept below 6 (plotext 6 is a full API rewrite and breaks the terminal plots).
- The `Procfile` runs `gunicorn vlbiplanobs.gui.main:server --workers 4 --timeout 120 ...`
  (see [Web Server](mode-server.md#production-deployment)).
- Package data (catalogs, templates, GUI assets) is installed correctly.

---

## Version 5.0.5

*Released 2026-10-02.*

### New Features

- **Network Monitoring Experiments** – `planobs ... --nme [--sched CODE]` plans NME schedules (all-antenna
  fringe-finder scans with periodic ftp fringe tests) and writes the NME `.key` file.
- **Observation reports** – `planobs observe -o FILE` saves all inputs and results as `.pdf`, `.txt`, `.md`, or `.json`.
- **SCHED templates** – `--template` selects a custom `.key` template and `planobs --get-key-template FILENAME`
  copies the bundled one; `--setup` sets the frequency setup line.
- Improved scheduler output for `--sched`.
- **GUI** – Redesigned observation-planning dashboard; copy button (Export to Polaris or copy link to clipboard);
  multi-source export to/from Polaris; fonts served locally (no third-party requests).
- CLI help is organized in groups and shown through a built-in pager.

### Improvements and Bug Fixes

- The former `--gui` and `--no-tui` options were removed; the observation output is always shown in the terminal.
- Durations are no longer rounded down to the nearest 10 minutes; short durations are allowed.
- Grouped antenna chips are properly sorted, and the antenna callbacks stay aligned with them.
- Fixed `borb` (PDF generation) installation on GNU/Linux and macOS.

## Version 5.0.4

*Released 2026-07-22.*

- Bug fix: removed a leftover `name` property in `sources.py`.

## Version 5.0.2

*Released 2026-07-22.* First 5.0 release. It includes all the changes listed under [Version 4.7.0](#version-470)
(grouped antenna chips, local horizons, terminal elevation plot, phase calibrator and check-source selection,
polarization calibration), plus:

- **Multiple sources** in the GUI, and export/import of the GUI configuration (including **Export to Polaris**).
- The GUI computes results in real time (no "compute" button).
- pySCHED station names are used in the schedules; faster `planobs source` visibility.
- 4 Gbps is now the default EVN data rate.
- GST ranges are shown again when an epoch is defined.
- File logging is opt-in via `--logging [LOGFILE]`.
- `kaleido` pinned below 1 (newer versions need Chrome to put figures in PDFs).
- Bug fixes: durations with more than one decimal fell back to 24 h; mandatory stations set to `all`; GUI networks
  after station changes; hardening against old cached values.

## Version 4.9.3

*Released 2026-04-08.*

- Python 3.12+ is now required.
- Handles observations where no baseline is visible.
- Fixed AstroGeo links for sources with negative declination.

## Version 4.9.2

*Released 2026-04-01.*

- Packaging fixes: missing CSS assets, `__version__`, log file path (`~` expansion), and relaxed `borb` version.

## Version 4.9

*Released 2026-03-30.*

- GUI layout changes and real-time dashboard mode.
- The terminal plot is limited to the targets, and unknown sources give a clearer error.

## Version 4.8.1

*Released 2026-03-16.*

### New Features

- **New CLI modes** – `planobs fringefinders`, `planobs phasecals` (calibrator searches in the RFC catalog),
  `planobs source` (source information, with `--gst` for GST observing windows), and `planobs antenna`
  (antenna information).
- First fully working version of the scheduler and `.key` file generation (`--sched`).
- Documentation pages for the new modes.

### Bug Fixes

- The Sun-constraint check was hard-coded to a given year.
- eMERLIN frequencies, and Ce/Ho antenna acceleration parameters.
- Dark mode and uv-plot highlighting in the GUI; short durations allowed in the GUI; style for small screens.

---

## Version 4.7.0

!!! note
    These changes were developed after 4.8.1 (under the 4.7.0 / 5.0a1 labels) and were first released in
    [Version 5.0.2](#version-502).

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

### From 5.0.x to 5.1.0

- Replace the renamed CLI options listed in the [5.1.0 migration table](#breaking-or-deprecated-cli-changes):
  `-t1`/`--starttime` → `-e`/`--epoch`, `--targets` → `-t`/`--target`, `fringefinders -t` → `-e`,
  `phasecals -n`/`--n-sources` → `-l`/`--max-lines`, `--catalog-file` → `--rfc-catalog`.
- Schedules generated with `--sched` may differ from 5.0.x, because added calibrators are now placed where they
  are visible.
- Python code that read `station.datarate` after building an observation should use
  `obs.station_datarate(codename)` or `obs.station_datarates` instead.
- `RFCCatalog.get_source()` now requires the exact name and returns `None` when the source is not found.

### From 4.5.x to 4.6.x

No breaking changes. New scheduling features are additive.

### From EVN Calculator

PlanObs is a complete rewrite. Key differences:

1. **Local Installation** – Can run offline
2. **Python API** – Full programmatic access
3. **Multiple Outputs** – GUI, CLI, and schedule files
4. **Extended Networks** – Support for non-EVN arrays
