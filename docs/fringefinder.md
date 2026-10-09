# Fringe Finder Sources

The **fringefinders** mode searches the RFC catalog for bright calibrator sources suitable for fringe detection during a VLBI observation. These sources are used to find and correct instrumental delays between telescopes.

## Usage

```bash
planobs fringefinders (-n NETWORK | -s STATIONS) -e EPOCH -d DURATION [OPTIONS]
```

### Quick Example

```bash
planobs fringefinders -s Ef Hh Mc Tr -e '2025-03-15 08:00' -d 8 -b 6cm
```

---

## Required Options

| Argument | Description |
|----------|-------------|
| `-n`, `--network` | One or more VLBI networks (e.g. `EVN`); their default stations are used. Can be combined with `-s`. |
| `-s`, `--stations` | List of antenna codenames or full names (e.g. `Ef Hh Mc Tr`). |
| `-e`, `--epoch` | Start of the observation in `'YYYY-MM-DD HH:MM'` format (UTC). |
| `-d`, `--duration` | Duration of the observation in hours. |

At least one of `-n` or `-s` is required.

!!! warning "Renamed option (v5.1.0)"
    The start time was previously given with `-t`/`--starttime`. `-t` is no longer accepted here (it means
    `--target` in the other subcommands) and produces an error; `--starttime` still works with a deprecation warning.

!!! note "Choosing fringe finders for a schedule"
    This mode only lists candidates. To put specific fringe finders in a schedule, use `--fringefinders` in the
    observe mode, which also accepts the `name/coordinates` format to override the catalog lookup:

    ```bash
    planobs -b 6cm -t 'M87' -n EVN -e '2025-03-15 08:00' -d 8 --sched eg123a \
        --fringefinders 'MyCal/12h30m49s +12d23m28s'
    ```

    The part before `/` is used as the name, the part after is parsed as coordinates.

---

## Optional Options

| Argument | Default | Description |
|----------|---------|-------------|
| `-b`, `--band` | all bands | Observing band (e.g. `6cm`, `18cm`). When provided, shows flux at that band only. |
| `--min-flux` | `0.5` | Minimum unresolved flux threshold in Jy. |
| `--min-elevation` | `20` | Minimum elevation in degrees. |
| `-l`, `--max-lines` | `20` | Maximum number of sources to return. |
| `--logging [LOGFILE]` | off | Log to a file. |
| `--require-all` | off | Require source to be visible by **all** stations, not just any. |
| `--station-catalog` | built-in | Path to a custom station catalog file. |
| `--json` | off | Output results in JSON format instead of a table. |

---

## Output

The tool prints a table with the following columns:

- **Name** – source name from the RFC catalog.
- **IVS Name** – International VLBI Service name.
- **Min elev. (deg)** – minimum elevation across all specified stations.
- **Total flux (Jy)** – total flux density at the observing band.
- **Unresolved (Jy)** – unresolved flux density (most relevant for fringe finding).
- **Bands** – bands where the source has been observed.
- **url** – link to the AstroGeo database page.
- **Antenna Visibility** – whether all or only some antennas see the source, and for all or part of the time.

### Visibility strip

Right below each source line there is a strip of coloured squares showing when the source can be observed
along the observation, from the start (left) to the end (right):

| Square | Meaning |
|--------|---------|
| Green | **All** antennas can observe the source. |
| Yellow | Not all antennas, but **more than 3**, can observe it. |
| Black | Any other case (3 or fewer antennas, and not all of them). |

```text
Antennas observing: ■ all  ■ more than 3  ■ fewer  (10 min each, from 08:00 UTC)

 Name                IVS Name   ...
 J1719+0817          1717+083   ...
 ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
```

- An antenna can observe the source when it is above `--min-elevation` and within the antenna limits
  (the same criterion used to select the candidates).
- The strip is as wide as your terminal allows. The time covered by one square is the shortest of
  1, 2, 5, 10, 15, 30 min, 1, 2, 3 or 6 h for which the whole observation fits; widen the terminal to get finer
  squares. The legend above the table shows the time per square and the start time.
- A square takes the colour of the worst moment within its time span: it is only green if all antennas see
  the source during the whole span.
- If the output has no colours (e.g. redirected to a file), the three levels are drawn as `■`, `□` and `·`.

When `--json` is used, a JSON object is printed with the search parameters, `total_found`, `shown`, and the
same data in the `sources` list (without the visibility strip).

---

## Examples

### EVN fringe finders at 6 cm

```bash
planobs fringefinders -s Ef Hh Mc Tr Wb O8 -e '2025-06-15 20:00' -d 12 -b 6cm
```

Finds candidates for a 12-hour EVN observation at 6 cm.

### Using a whole network

```bash
planobs fringefinders -n EVN -e '2025-06-15 20:00' -d 6.5 -b 18cm
```

### High elevation and flux requirements

```bash
planobs fringefinders -s Ef Hh Mc Tr -e '2025-06-15 20:00' -d 8 -b 6cm \
    --min-elevation 30 --min-flux 0.8
```

Returns only sources above 30° elevation with ≥ 0.8 Jy unresolved flux.

### Strict: visible by all stations, top 5 only

```bash
planobs fringefinders -s Ef Hh Mc Tr Ys -e '2025-06-15 20:00' -d 8 -b 6cm \
    --require-all --min-flux 1.0 -l 5
```

### JSON output for scripting

```bash
planobs fringefinders -s Ef Hh Mc Tr -e '2025-03-15 08:00' -d 8 -b 6cm --json
```

---

## Tips for Good Fringe Finders

1. **Flux** – look for sources with high unresolved flux (> 0.5 Jy is typical).
2. **Elevation** – higher elevation means less atmospheric absorption.
3. **Visibility** – the source should be above the horizon for all stations during the observation.
4. **Sky distribution** – sources spread across the sky improve delay calibration.
5. **Compactness** – point-like sources (high unresolved/total ratio) are ideal.

---

## Troubleshooting

**No sources found** – lower `--min-flux`, reduce `--min-elevation`, or remove `--require-all`.

**Too many sources** – increase `--min-flux`, add `--require-all`, or use `-l` to limit output.

**Station not found** – check the codename with `planobs --list-antennas`.
