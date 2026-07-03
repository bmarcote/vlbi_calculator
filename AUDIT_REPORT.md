# Audit Report: planobs (whole program)

Scope: full `src/vlbiplanobs` package (GUI, CLI, computation core)
Date: 2026-07-03 (updated same day: all remaining findings fixed)
Files in scope: 18 Python modules (~12,000 lines)

## Summary

The program is in good shape overall: the computation core is consistent, the GUI callbacks are
well-factored, and the 93-test suite passes. One production crash (the Dash callback `TypeError`
reported from the server logs), nine latent correctness bugs, and six low-severity maintenance
items were found and **all fixed in this audit**; every fix is verified by the test suite plus
targeted end-to-end checks. No security issues were found. No open findings remain.

## Fixed during this audit

### 1. Dash callback crash on every menu-highlight dispatch (the reported server error)

**Location**: `src/vlbiplanobs/gui/callbacks.py:225`
**Category**: Logic
**Problem**: `highlight_active_menu_item` declares a single (non-list) `Output` with an `ALL`
pattern, so `ctx.outputs_list` *is* the flat list of output specs. The code iterated
`ctx.outputs_list[0]` — the first spec dict — yielding its string keys and crashing with
`TypeError: string indices must be integers, not 'str'` on every dispatch with group chips present.
**Fix**: iterate `ctx.outputs_list` directly (the two other pattern-matching callbacks declare
list-wrapped outputs/states, so their `[0]` indexing is correct).

### 2. Station mount axis parameters kept as raw config strings

**Location**: `src/vlbiplanobs/stations.py:982-990`
**Category**: Logic
**Problem**: `_parse_station_from_configfile` passed `station["ax2acc"]` (and, in the no-acceleration
branch, `station["ax1rate"]`/`station["ax2rate"]`) — raw strings from the catalog — into `Axis`
instead of the parsed `configs[...]` Quantities. Confirmed live: Effelsberg's axis-2 acceleration
was the string `'0.018'`. Any code using `Axis.speed`/`Axis.acceleration` arithmetic
(`Station.slewing_time`) would raise `TypeError`.
**Fix**: use the parsed `configs` values in both branches.

### 3. Filtered station sets corrupt per-band max data rates

**Location**: `src/vlbiplanobs/stations.py:766, 812`
**Category**: Logic
**Problem**: `filter_antennas` and `filter_networks` passed `self._bands.values()` (a `dict_values`
object) as `max_datarates`. `Stations.__init__` only treats a `Sequence` as per-band values, so the
`dict_values` object was assigned wholesale as the max datarate of *every* band. Confirmed live:
`filter_antennas([...]).max_datarate('18cm')` returned `dict_values([None, ...])` instead of `None`.
**Fix**: wrap in `list(...)`.

### 4. `Observation.sources(source_type=...)` filter poisons/ignores the cache

**Location**: `src/vlbiplanobs/observation.py:290-316`
**Category**: Logic
**Problem**: the full source list was cached on first call regardless of the `source_type` argument:
a first filtered call cached a filtered list served to later unfiltered callers, and a first
unfiltered call made later filtered calls return unfiltered results.
**Fix**: always cache the full list; apply the filter on return.

### 5. `is_always_observable()` cache never invalidated

**Location**: `src/vlbiplanobs/observation.py` (scans/times/duration/stations setters)
**Category**: Logic / Consistency
**Problem**: the setters cleared a phantom attribute `_is_always_observable`, while the actual cache
attribute is `_is_always_visible` — so changing scans, times, duration, or stations on an existing
`Observation` returned stale always-observable results. The `_source_list` cache was similarly never
cleared when scans changed.
**Fix**: setters now clear `_is_always_visible`, and the scans setter also clears `_source_list`.

### 6. `baseline_sensitivity()` bakes the first integration time into the cache

**Location**: `src/vlbiplanobs/observation.py:1290-1340`
**Category**: Logic
**Problem**: the per-baseline sensitivities were divided by `sqrt(integration_time)` *in place* on
first call, so any later call with a different `integration_time` silently returned values scaled
for the first one. Additionally, the single-antenna branches looked up keys in both orientations
(`A-B` and `B-A`) although only one orientation is stored — a guaranteed `KeyError`.
**Fix**: cache raw (1-second) values, scale on every return; single-antenna lookup now matches keys
by codename membership.

### 7. GST range end minutes printed from the start time

**Location**: `src/vlbiplanobs/observation.py:1640-1643` (`print_obs_times`)
**Category**: Logic
**Problem**: the end-of-range minutes used `gstimes[0]` instead of `gstimes[-1]`, so the printed GST
range always ended with the start time's minutes (e.g. "06:44-14:44" instead of "06:44-14:45").
**Fix**: use `gstimes[-1]` for the end minutes.

### 8. Shared cache-year between the two module reference-time caches

**Location**: `src/vlbiplanobs/observation.py:26-47, 153-172`
**Category**: Logic
**Problem**: `_REF_TIMES` and `_REF_YEAR` shared a single `_CACHE_YEAR`. On a year rollover in a
long-running server, whichever property was hit first bumped the year and the other cache was never
refreshed, serving reference times for the previous year indefinitely.
**Fix**: each cache now tracks its own year (`_CACHE_YEAR_TIMES` / `_CACHE_YEAR_YEAR`).

### 9. `ScanBlock.fractional_time()` wrong with multiple distinct `every` values

**Location**: `src/vlbiplanobs/sources.py:854-868`
**Category**: Logic
**Problem**: over `lcm(every_i)` cycles, an `every=N` scan is observed `lcm/N` times, but the total
duration summed `duration*(N/lcm)` — correct only when a single distinct `every` value exists.
With e.g. `every=2` and `every=4` scans, fractions did not sum to 1. Also `every == 0` caused a
`ZeroDivisionError`.
**Fix**: total now uses `duration*(lcm/N)`; the divisor guard uses `every > 0`. Verified that
fractions sum to 1.0 for mixed `every` values.

### 10. `ant_warning` only reports the last source's excluded antennas

**Location**: `src/vlbiplanobs/gui/outputs.py:243-245`
**Category**: Logic
**Problem**: the loop overwrote `ants_excluded` on each source instead of accumulating, so for
multi-source observations only the last source's non-observing antennas were shown.
**Fix**: accumulate (deduplicated) across sources.

## Low-severity findings (fixed in the follow-up pass, same day)

### L1. Temporary files from PDF export were never deleted — FIXED

**Location**: `src/vlbiplanobs/gui/outputs.py` (`summary_pdf`), `src/vlbiplanobs/gui/main.py`
(`download_pdf_per_target`)
**Category**: Best practices
**Problem**: `NamedTemporaryFile(delete=False)` PNG/PDF files accumulated in the server's temp dir.
**Fix**: `summary_pdf` now removes the figure PNG after the PDF is written (the PNG path is
captured before rendering so it is cleaned up even if figure export fails); the download callback
removes the PDF after `dcc.send_file` has read it into the response. Verified: no PNGs remain
after generation.

### L2. `url_open` crashed the callback on malformed `config` URL parameter — FIXED

**Location**: `src/vlbiplanobs/gui/callbacks.py` (`url_open`)
**Category**: Logic (robustness)
**Problem**: a hand-edited or truncated `?config=` query string raised in `json.loads`/`Version`,
producing a 500 on that dispatch.
**Fix**: version/config parsing wrapped in try/except raising `PreventUpdate`; verified with
malformed JSON, empty config list, and an invalid version string, plus a valid config still parsing.

### L3. `functools.cache` on `SourceCatalog` instance methods — FIXED

**Location**: `src/vlbiplanobs/sources.py` (`SourceCatalog.source_names`/`.sources`)
**Category**: Best practices
**Problem**: caches keyed on `self` kept catalogs alive (memory) and went stale if
`read_personal_catalog` was called after first use.
**Fix**: replaced with per-instance dict caches (keyed by `include_calibrators`) that
`read_personal_catalog` clears before re-reading. Verified against `tests/targets.toml`.

### L4. Multi-block elevation heatmap used the last block's y-axis labels — FIXED

**Location**: `src/vlbiplanobs/gui/plots.py` (`elevation_plot`, `elevation_plot_from_data`)
**Category**: Logic (latent)
**Problem**: `n_ants`/`ant_names` from the last loop iteration configured only the first subplot's
y-axis (other subplots got default numeric ticks).
**Fix**: axis styling is now applied to all subplots via `update_yaxes`, and each subplot gets its
own block's antenna names as tick labels. Verified with a two-target observation: both `yaxis` and
`yaxis2` carry the correct labels.

### L5. Inconsistent `--min-flux` default in `main_phasecal` — FIXED

**Location**: `src/vlbiplanobs/calibrators.py` (`main_phasecal`), `src/vlbiplanobs/cli.py`
(`handle_phase_cal_command`)
**Category**: Consistency
**Problem**: `main_phasecal` defaulted to 0.0 while its help text and the `planobs phasecals`
wrapper documented 0.1; the wrapper only forwarded non-default values, so the effective default
silently diverged from the docs.
**Fix**: `main_phasecal` now defaults to the documented 0.1, and the wrapper always forwards
`--min-flux` so the two parsers cannot drift again.

### L6. `VLBIObs.summary(gui=False, tui=False)` silently returned None — FIXED

**Location**: `src/vlbiplanobs/cli_obs.py` (`VLBIObs.summary`)
**Category**: Consistency
**Fix**: explicit `return None` and a docstring noting the behavior (the CLI warns about this
combination before calling).

## No findings

- **Security**: no injection risks found. The web GUI takes only typed/validated inputs (source
  specs are validated through `astropy`/catalog lookups before use; upload contents are
  base64-decoded with explicit error handling and never executed or written to disk; the SCHED
  key-file writer only runs locally via the CLI). No secrets are handled or logged.
- **Concurrency**: the `Observation` caches are guarded by `RLock`s; the per-target
  `ThreadPoolExecutor` in the GUI builds independent observation objects.
- **Dependency hygiene**: pinned via `uv.lock`; no known-vulnerable patterns spotted in use.

## Verification

- Full test suite: 93/93 passed after all fixes (including the L1–L6 follow-up pass).
- Ruff (`E9,F`): clean.
- GUI app import + callback registration: OK.
- CLI end-to-end (`planobs -b 6cm -n EVN -t 3C273 -d 4`, `planobs phasecals -t J1229+0203 -b 6cm`): OK.
- Simulated Dash dispatch of the fixed callback: returns correct class list.
- PDF generation: no leaked temp PNGs; PDF removed after serving.
- `url_open`: malformed config/version parameters ignored; valid configs still applied.
- Two-target elevation heatmap: both subplots correctly labeled.
