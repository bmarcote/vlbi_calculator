# CLI Reference

Complete reference for the `planobs` command-line interface.

## Basic Usage

```bash
planobs -b BAND [OPTIONS]
```

## Required Arguments

| Argument | Description |
|----------|-------------|
| `-b`, `--band` | Observing band (e.g., `6cm`, `18cm`, `1.3cm`) |

## Source Options

| Argument | Description |
|----------|-------------|
| `-t`, `--targets` | Source name(s) or coordinates. Resolves via SIMBAD/NED |
| `-sc`, `--source-catalog` | Path to custom source catalog file |

## Station Options

| Argument | Description |
|----------|-------------|
| `-n`, `--network` | VLBI network(s): `EVN`, `eMERLIN`, `VLBA`, `LBA`, etc. |
| `-s`, `--stations` | Individual station codes (e.g., `Ef`, `Wb`, `Jb`) |
| `--station-catalog` | Path to custom station catalog file |

## Time Options

| Argument | Description |
|----------|-------------|
| `-t1`, `--starttime` | Observation start: `'YYYY-MM-DD HH:MM'` (UTC) |
| `-d`, `--duration` | Observation duration in hours |

## Scheduling Options

| Argument | Description |
|----------|-------------|
| `--sched` | Generate pySCHED `.key` file with given name |
| `--fringefinders` | Fringe finder sources or count (default: `2`) |
| `--polcal` | Include polarization calibration |
| `--pulsar` | Include pulsar source |

## Output Options

| Argument | Description |
|----------|-------------|
| `--data-rate` | Maximum data rate in Mb/s |
| `--gui` | Show plots in browser |
| `--no-tui` | Suppress terminal output |
| `--debug` | Show debug messages |

## Information Commands

| Argument | Description |
|----------|-------------|
| `--list-antennas` | List all available antennas |
| `--list-networks` | List all VLBI networks |
| `--list-bands` | List all observing bands |

## Examples

### Basic Visibility Check

```bash
planobs -b 6cm -t 'M87' --network EVN
```

### Full Observation Planning

```bash
planobs -b 6cm \
  -t 'Cygnus A' \
  --network EVN eMERLIN \
  --starttime '2025-06-15 08:00' \
  --duration 12 \
  --data-rate 2048
```

### Generate Schedule File

```bash
planobs -b 6cm \
  -t 'M87' \
  --network EVN \
  --starttime '2025-03-15 08:00' \
  --duration 8 \
  --sched eg123a
```

### Use Custom Catalogs

```bash
planobs -b 6cm \
  --source-catalog my_sources.inp \
  --station-catalog my_stations.inp \
  -t 'MyTarget'
```

### List Available Resources

```bash
planobs --list-networks
planobs --list-antennas
planobs --list-bands
```
