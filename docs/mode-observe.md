# Observation Planning

The **observe** mode is the core of PlanObs. It plans a VLBI observation by computing source visibility, estimating sensitivity, and optionally generating a pySCHED-compatible schedule file.

## Usage

```bash
planobs observe -b BAND [OPTIONS]
```

!!! tip "Legacy syntax"
    You can omit the `observe` subcommand. Running `planobs -b 6cm ...` is equivalent to `planobs observe -b 6cm ...`.

---

## Required Arguments

| Argument | Description |
|----------|-------------|
| `-b`, `--band` | Observing band (e.g. `6cm`, `18cm`, `1.3cm`). Use `--list-bands` to see all options. |

You must also specify at least one of `--network` or `--stations`.

---

## Station Selection

| Argument | Description |
|----------|-------------|
| `-n`, `--network` | One or more VLBI network names (e.g. `EVN`, `eMERLIN`, `VLBA`, `LBA`). Takes the default stations of each network. |
| `-s`, `--stations` | Individual station codenames or full names (e.g. `Ef`, `Wb`, `Jb`). Can be combined with `--network`. |
| `--station-catalog` | Path to a custom station catalog file. Overrides the built-in catalog. |

### Examples

Use a predefined network:

```bash
planobs -b 6cm --network EVN
```

Combine a network with extra stations:

```bash
planobs -b 6cm --network EVN --stations Ar Ys
```

Use only specific stations:

```bash
planobs -b 6cm --stations Ef Hh Mc Tr Wb On
```

---

## Source Selection

| Argument | Description |
|----------|-------------|
| `-t`, `--targets` | One or more source names or coordinates. Names are resolved via SIMBAD/NED/VizieR. Coordinates use `hh:mm:ss dd:mm:ss` or `XXhXXmXXs XXdXXmXXs` format. |
| `-sc`, `--source-catalog` | Path to a TOML source catalog file. When combined with `--targets`, only the named blocks from the file are used. |

### Examples

Observe a named source:

```bash
planobs -b 6cm -t 'M87' --network EVN
```

Observe a source by coordinates:

```bash
planobs -b 6cm -t '12h30m49.4s +12d23m28s' --network EVN
```

Use a source catalog:

```bash
planobs -b 6cm --source-catalog my_sources.toml --network EVN
```

Select specific blocks from the catalog:

```bash
planobs -b 6cm --source-catalog my_sources.toml -t block1 block2 --network EVN
```

---

## Time and Duration

| Argument | Description |
|----------|-------------|
| `-t1`, `--starttime` | Start of the observation in `'YYYY-MM-DD HH:MM'` format (UTC). |
| `-d`, `--duration` | Total duration of the observation in hours. |

When **no time is specified**, PlanObs searches for the optimal GST range and reports when the source is observable.

When **both start time and duration are given**, PlanObs evaluates visibility and sensitivity for that specific epoch.

!!! warning
    If you provide `--starttime`, you must also provide `--duration`.

### Examples

No specific time (find optimal GST range):

```bash
planobs -b 6cm -t 'M87' --network EVN
```

Fixed epoch observation:

```bash
planobs -b 6cm -t 'M87' --network EVN --starttime '2025-06-15 08:00' --duration 12
```

---

## Setup Options

| Argument | Description |
|----------|-------------|
| `--data-rate` | Maximum data rate in Mb/s. If omitted, PlanObs infers it from the network. |

---

## Scheduling

| Argument | Description |
|----------|-------------|
| `--sched` | Generate a pySCHED `.key` schedule file. The value is used as the experiment code and filename. |
| `--fringefinders` | Fringe finder source(s) or a number of automatic selections (default: `2`). |
| `--polcal` | Include polarization calibration scans. |
| `--pulsar` | Include pulsar source scans. Accepts a source name or a number referencing the source catalog. |

### Example

Generate a schedule file with automatic fringe finder selection:

```bash
planobs -b 6cm -t 'M87' --network EVN \
  --starttime '2025-03-15 08:00' --duration 8 \
  --sched eg123a
```

This creates `eg123a.key`. See **[Scheduling](scheduling.md)** for details on the `.key` file format.

---

## Output Options

| Argument | Description |
|----------|-------------|
| `--gui` | Open graphical plots in the browser. |
| `--no-tui` | Suppress the terminal (TUI) output. |
| `--debug` | Show debug messages and execution time. |

---

## Information Commands

These print reference data and exit immediately:

| Argument | Description |
|----------|-------------|
| `--list-antennas` | List all antennas in the catalog with their bands and locations. |
| `--list-networks` | List all known VLBI networks and their default stations. |
| `--list-bands` | List all observing bands and which networks support them. |

```bash
planobs --list-networks
planobs --list-antennas
planobs --list-bands
```

---

## Output

When run with a target source, PlanObs prints:

1. **Observation summary** – band, duration, stations, setup (data rate, bandwidth, subbands).
2. **Source list** – all sources in every scan block, with coordinates.
3. **Visibility chart** – a text-based chart showing when each station can observe the source.
4. **Optimal observing window** – GST or UTC range when the maximum number of stations can observe.
5. **Thermal noise estimate** – expected RMS noise for the target, with on-source time.
6. **Sun constraint warning** – if the Sun is too close to the source during the observation.

When run without a target (only band + stations + time), PlanObs reports the expected thermal noise for a hypothetical source at ±45° elevation.

---

## Full Examples

### Quick visibility check (no fixed time)

```bash
planobs -b 6cm -t 'Cygnus A' --network EVN eMERLIN
```

### Detailed observation at a fixed epoch

```bash
planobs -b 6cm \
  -t 'Cygnus A' \
  --network EVN eMERLIN \
  --stations Ar Ys \
  --starttime '2025-06-15 08:00' \
  --duration 12 \
  --data-rate 2048
```

### Generate a full schedule

```bash
planobs -b 6cm \
  -t 'M87' \
  --network EVN \
  --starttime '2025-03-15 08:00' \
  --duration 8 \
  --sched eg123a \
  --fringefinders 3 \
  --polcal
```

### Custom catalogs

```bash
planobs -b 6cm \
  --source-catalog my_sources.toml \
  --station-catalog my_stations.inp \
  -t MyTarget \
  --network EVN
```

### Terminal-only output (no GUI)

```bash
planobs -b 6cm -t 'M87' --network EVN --no-tui
```
