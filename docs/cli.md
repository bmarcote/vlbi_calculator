# CLI Reference

The `planobs` command-line interface is organised into **modes** (subcommands). Each mode addresses a specific task in the VLBI observation planning workflow.

## Modes at a Glance

| Command | Purpose | Key arguments |
|---------|---------|---------------|
| `planobs [observe]` | Plan a VLBI observation | `-b BAND`, `-t TARGET`, `-n NETWORK` |
| `planobs fringefinders` | Find fringe finder sources | `-n NETWORK`/`-s STATIONS`, `-e EPOCH`, `-d DURATION` |
| `planobs phasecals` | Find phase calibrator sources | `TARGET` |
| `planobs source` | Look up source information | `<source_name>`, `--gst` |
| `planobs antenna` (alias `ant`) | Look up antenna information | `[ANTENNA]`, `-b BAND` |
| `planobs server` | Launch the web GUI | `--host`, `--port` |

!!! tip "Legacy syntax"
    Running `planobs -b 6cm ...` (without a subcommand) is equivalent to `planobs observe -b 6cm ...`.

---

## Quick Examples

### Plan an observation

```bash
planobs -b 6cm -t 'M87' --network EVN
```

### Find fringe finders

```bash
planobs fringefinders -s Ef Hh Mc Tr -e '2025-03-15 08:00' -d 8 -b 6cm
```

### Find phase calibrators

```bash
planobs phasecals 'M87' -b 6cm
```

### Look up a source

```bash
planobs source '3C273'
```

### Look up an antenna

```bash
planobs antenna Ef
planobs antenna -b 1.3cm
```

### Start the web server

```bash
planobs server
```

---

## Option Names

The same concept uses the same flag in every subcommand:

| Option | Meaning | Used in |
|--------|---------|---------|
| `-b`, `--band` | Observing band | observe, fringefinders, phasecals, antenna |
| `-n`, `--network` | VLBI network(s) | observe, fringefinders |
| `-s`, `--stations` | Individual stations | observe, fringefinders |
| `-e`, `--epoch` | Start time, `'YYYY-MM-DD HH:MM'` (UTC) | observe, fringefinders |
| `-d`, `--duration` | Duration in hours (float) | observe, fringefinders |
| `-t`, `--target` | Target source(s) | observe, phasecals |
| `-l`, `--max-lines` | Maximum number of results | fringefinders, phasecals |
| `-sc`, `--source-catalog` | Personal source catalog | observe, phasecals |
| `--station-catalog` | Personal station catalog | observe, fringefinders |
| `--json` | Machine-readable output | fringefinders, phasecals |
| `--logging [LOGFILE]` | Log to a file | all |

!!! warning "Renamed options (v5.1.0)"
    Deprecated spellings still work but print a warning: `-t1`/`--starttime` (use `-e`/`--epoch`),
    `--targets` (use `-t`/`--target`), `--n-sources` (use `-l`/`--max-lines`) and `--catalog-file`
    (use `--rfc-catalog`). Two short options were removed because they now mean something else:
    `planobs fringefinders -t` (use `-e`) and `planobs phasecals -n` (use `-l`).

---

## Detailed Mode References

Each mode has its own dedicated documentation page with full argument tables, output descriptions, and worked examples:

- **[Observation Planning](mode-observe.md)** – `planobs [observe]`
- **[Fringe Finders](fringefinder.md)** – `planobs fringefinders`
- **[Phase Calibrators](phasecal.md)** – `planobs phasecals`
- **[Source Lookup](mode-source.md)** – `planobs source` (also covers `planobs antenna`, see [Antenna Lookup](mode-source.md#antenna-lookup-planobs-antenna))
- **[Web Server](mode-server.md)** – `planobs server`

---

## Common Patterns

### List available resources

These flags work in the default (observe) mode and print reference data:

```bash
planobs --list-networks   # all known VLBI networks
planobs --list-antennas   # all antennas with bands and locations
planobs --list-bands      # all observing bands and supporting networks
```

### Custom catalogs

The observe mode accepts custom source and station catalogs (fringefinders accepts `--station-catalog`, phasecals accepts `-sc` and `--rfc-catalog`):

```bash
planobs -b 6cm --source-catalog my_sources.toml --station-catalog my_stations.inp \
  -t MyTarget --network EVN
```

### Observation reports

The observe mode accepts `-o`/`--output` with a `.pdf`, `.txt`, `.md`, or `.json` filename. The extension selects the format; JSON includes both reusable inputs and calculated outputs.

```bash
planobs observe -b 6cm -t 'M87' --network EVN -o m87-summary.pdf
planobs observe -b 6cm -t 'M87' --network EVN -o m87-results.json
```

The fringefinders and phasecals modes separately support `--json` for machine-readable terminal output:

```bash
planobs fringefinders -s Ef Hh Mc Tr -e '2025-03-15 08:00' -d 8 -b 6cm --json
planobs phasecals 'M87' -b 6cm --json
```

### Generate a schedule file

```bash
planobs -b 6cm -t 'M87' --network EVN \
  --epoch '2025-03-15 08:00' --duration 8 --sched eg123a
```

### Generate a schedule file with a given frequency setup

```bash
planobs -b 6cm -t 'M87' --network EVN \
  --epoch '2025-03-15 08:00' --duration 8 --sched eg123a \
  --setup 'evn6cm-2Gbps-32MHz.set'
```

### Auto-select calibrators in schedule

```bash
planobs -b 6cm -t 'M87' --network EVN \
  --epoch '2025-03-15 08:00' --duration 8 --sched eg123a \
  --fringefinders 3 --phasecal --check-source
```

See **[Scheduling](scheduling.md)** for details on the `.key` file format and auto-selection features.
