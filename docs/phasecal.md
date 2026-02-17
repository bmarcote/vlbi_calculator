# Phase Calibrator Sources

The **phasecals** mode searches the RFC catalog for compact calibrator sources near your target, suitable for phase-referenced VLBI observations. Phase calibrators correct atmospheric and instrumental phase variations.

## Usage

```bash
planobs phasecals -t TARGET [OPTIONS]
```

### Quick Example

```bash
planobs phasecals -t 'M87' -b 6cm
```

---

## Required Options

| Argument | Description |
|----------|-------------|
| `-t`, `--target` | Target source name. Accepts J2000 names, IVS names from the RFC catalog, or any name resolvable by SIMBAD/NED/VizieR. |

---

## Optional Options

| Argument | Default | Description |
|----------|---------|-------------|
| `-b`, `--band` | all bands | Observing band (e.g. `6cm`, `1.3cm`). When omitted, shows average flux across all bands. |
| `--max-separation` | `5.0` | Maximum angular separation from the target in degrees. |
| `--min-flux` | `0.1` | Minimum unresolved flux threshold in Jy. |
| `-n`, `--n-sources` | all | Maximum number of sources to return (sorted by separation). |
| `--catalog-file` | built-in | Path to a custom RFC catalog file. |
| `--json` | off | Output results in JSON format instead of a table. |

---

## Output

The tool prints a table with the following columns:

- **Name** – source name from the RFC catalog.
- **IVS Name** – International VLBI Service name.
- **Separation (deg)** – angular distance from the target.
- **Total Flux (Jy)** – total flux density at the observing band.
- **Unresolved (Jy)** – unresolved flux density (most relevant for phase calibration).
- **Bands** – bands where the source has been observed.
- **url** – link to the AstroGeo database page.

When `--json` is used, the same data is printed as a JSON array.

---

## Examples

### Basic search at 6 cm

```bash
planobs phasecals -t 'M87' -b 6cm
```

Finds all calibrator candidates within 5° of M87 at 6 cm.

### Tight separation constraint

```bash
planobs phasecals -t 'M87' -b 6cm --max-separation 2 --min-flux 0.2
```

Calibrators within 2° with ≥ 0.2 Jy unresolved flux.

### All-band search, top 10

```bash
planobs phasecals -t '3C273' --max-separation 10 -n 10
```

Searches all bands for calibrators within 10° of 3C273, returns the 10 closest.

### High-frequency observation

```bash
planobs phasecals -t 'NGC 4258' -b 1.3cm --max-separation 3 --min-flux 0.5
```

Stricter flux requirements for a 1.3 cm observation.

### JSON output for scripting

```bash
planobs phasecals -t 'M87' -b 6cm --json
```

---

## Phase Referencing Guidelines

### Separation distance

| Separation | Quality |
|------------|---------|
| < 1° | Excellent — minimal phase decorrelation |
| 1–3° | Good — commonly used range |
| 3–5° | Acceptable at lower frequencies |
| > 5° | Risk of significant decorrelation, especially at high frequencies |

### Flux requirements by frequency

| Frequency range | Recommended minimum flux |
|-----------------|--------------------------|
| ≥ 6 cm | 0.1–0.2 Jy |
| 3–6 cm | 0.2–0.5 Jy |
| ≤ 2 cm | 0.5–1.0+ Jy |

### Source compactness

A high **unresolved/total flux ratio** indicates a compact source. Extended structure introduces phase errors — always check the AstroGeo link for source maps.

---

## Tips

1. **Proximity** – choose the closest suitable calibrator.
2. **Flux** – ensure sufficient flux for your band and integration time.
3. **Compactness** – prefer high unresolved/total flux ratios.
4. **Elevation** – consider the calibrator elevation during your observation.
5. **Backup** – identify at least one backup calibrator.

---

## Troubleshooting

**No sources found** – increase `--max-separation`, lower `--min-flux`, or verify the target name.

**Target not found** – check spelling, try coordinates, or verify with `planobs source <name>`.

**Too many weak sources** – increase `--min-flux`, reduce `--max-separation`, or use `-n` to limit output.

---

## Integration with Observation Planning

Once you have identified calibrators, add them to your source catalog file and include them in the observation plan:

```bash
planobs -b 6cm --source-catalog my_sources.toml --network EVN \
  --starttime '2025-06-15 08:00' --duration 8
```

PlanObs automatically computes phase-referencing cycle times based on the observing band.
