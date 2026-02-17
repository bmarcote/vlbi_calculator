# Fringe Finder Sources

The **fringefinders** mode searches the RFC catalog for bright calibrator sources suitable for fringe detection during a VLBI observation. These sources are used to find and correct instrumental delays between telescopes.

## Usage

```bash
planobs fringefinders -s STATIONS -t STARTTIME -d DURATION [OPTIONS]
```

### Quick Example

```bash
planobs fringefinders -s Ef Hh Mc Tr -t '2025-03-15 08:00' -d 8 -b 6cm
```

---

## Required Options

| Argument | Description |
|----------|-------------|
| `-s`, `--stations` | List of antenna codenames or full names (e.g. `Ef Hh Mc Tr`). |
| `-t`, `--starttime` | Start of the observation in `'YYYY-MM-DD HH:MM'` format (UTC). |
| `-d`, `--duration` | Duration of the observation in hours. |

---

## Optional Options

| Argument | Default | Description |
|----------|---------|-------------|
| `-b`, `--band` | all bands | Observing band (e.g. `6cm`, `18cm`). When provided, shows flux at that band only. |
| `--min-flux` | `0.5` | Minimum unresolved flux threshold in Jy. |
| `--min-elevation` | `20` | Minimum elevation in degrees. |
| `-l`, `--max-lines` | `20` | Maximum number of sources to return. |
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

When `--json` is used, the same data is printed as a JSON array.

---

## Examples

### EVN fringe finders at 6 cm

```bash
planobs fringefinders -s Ef Hh Mc Tr Wb On -t '2025-06-15 20:00' -d 12 -b 6cm
```

Finds candidates for a 12-hour EVN observation at 6 cm.

### High elevation and flux requirements

```bash
planobs fringefinders -s Ef Hh Mc Tr -t '2025-06-15 20:00' -d 8 -b 6cm \
    --min-elevation 30 --min-flux 0.8
```

Returns only sources above 30° elevation with ≥ 0.8 Jy unresolved flux.

### Strict: visible by all stations, top 5 only

```bash
planobs fringefinders -s Ef Hh Mc Tr Ys -t '2025-06-15 20:00' -d 8 -b 6cm \
    --require-all --min-flux 1.0 -l 5
```

### JSON output for scripting

```bash
planobs fringefinders -s Ef Hh Mc Tr -t '2025-03-15 08:00' -d 8 -b 6cm --json
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
