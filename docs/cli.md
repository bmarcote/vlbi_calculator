# CLI Reference

The `planobs` command-line interface is organised into **modes** (subcommands). Each mode addresses a specific task in the VLBI observation planning workflow.

## Modes at a Glance

| Command | Purpose | Key arguments |
|---------|---------|---------------|
| `planobs [observe]` | Plan a VLBI observation | `-b BAND`, `-t TARGET`, `-n NETWORK` |
| `planobs fringefinders` | Find fringe finder sources | `-s STATIONS`, `-t STARTTIME`, `-d DURATION` |
| `planobs phasecals` | Find phase calibrator sources | `-t TARGET` |
| `planobs source` | Look up source information | `<source_name>` |
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
planobs fringefinders -s Ef Hh Mc Tr -t '2025-03-15 08:00' -d 8 -b 6cm
```

### Find phase calibrators

```bash
planobs phasecals -t 'M87' -b 6cm
```

### Look up a source

```bash
planobs source '3C273'
```

### Start the web server

```bash
planobs server
```

---

## Detailed Mode References

Each mode has its own dedicated documentation page with full argument tables, output descriptions, and worked examples:

- **[Observation Planning](mode-observe.md)** – `planobs [observe]`
- **[Fringe Finders](fringefinder.md)** – `planobs fringefinders`
- **[Phase Calibrators](phasecal.md)** – `planobs phasecals`
- **[Source Lookup](mode-source.md)** – `planobs source`
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

Both the observe and fringefinders modes accept custom catalogs:

```bash
planobs -b 6cm --source-catalog my_sources.toml --station-catalog my_stations.inp \
  -t MyTarget --network EVN
```

### JSON output

The fringefinders and phasecals modes support `--json` for machine-readable output:

```bash
planobs fringefinders -s Ef Hh Mc Tr -t '2025-03-15 08:00' -d 8 -b 6cm --json
planobs phasecals -t 'M87' -b 6cm --json
```

### Generate a schedule file

```bash
planobs -b 6cm -t 'M87' --network EVN \
  --starttime '2025-03-15 08:00' --duration 8 --sched eg123a
```

See **[Scheduling](scheduling.md)** for details on the `.key` file format.
