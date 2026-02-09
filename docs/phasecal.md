# Phase Calibrator Sources

The `planobs_phasecal` tool helps you find phase calibrator sources near your target for phase-referenced VLBI observations. Phase calibrators are bright, compact sources used to correct atmospheric and instrumental phase variations.

## Overview

Phase calibrators must be:
- Close to your target source (typically <10 degrees)
- Bright enough to detect on short timescales
- Compact (point-like) to avoid structure phase errors
- Observable at the same time as your target

## Usage

### Basic Usage

```bash
planobs_phasecal -t 'M87' -b 6cm
```

### Advanced Usage

```bash
planobs_phasecal -t 'M87' -b 6cm --max-separation 3 --min-flux 0.2 -n 5
```

## Command Line Options

### Required Options

- `-t, --target`: Target source name
  - Can be a J2000 coordinate name from the RFC catalog
  - Can be an IVS name from the RFC catalog
  - Can be coordinates in various formats (handled by the main planobs tool)

### Optional Options

- `-b, --band`: Observing band (optional)
  - Supported bands: '18cm', '21cm', '13cm', '6cm', '5cm', '3.6cm', '2cm', '1.3cm', '0.7cm'
  - If not provided, the tool searches all bands and shows average fluxes
  - When provided, shows flux information specific to that band

- `--max-separation`: Maximum angular separation in degrees (default: 5.0)
  - Searches for calibrators within this distance from the target
  - Smaller values find closer calibrators (better for phase calibration)
  - Larger values provide more options but may be less optimal

- `--min-flux`: Minimum unresolved flux threshold in Jy (default: 0.1)
  - Higher values return fewer, brighter sources
  - Lower values return more candidates but may include weaker sources
  - For high-frequency observations, you may want higher flux thresholds

- `-n, --n-sources`: Maximum number of sources to return (default: all)
  - Limits output to the top N candidates by separation
  - Useful when you want to compare a few options

- `--catalog-file`: Path to custom RFC catalog file
  - Use if you have a modified or updated RFC catalog
  - Overrides the default catalog included with PlanObs

## Output

The tool displays a table with the following columns:

- **Name**: Source name from the RFC catalog
- **IVS Name**: International VLBI Service name
- **Separation (deg)**: Angular separation from the target source
- **Total Flux (Jy)**: Total flux density at the observing band
- **Unresolved (Jy)**: Unresolved flux density (most relevant for phase calibration)
- **Bands**: Bands where the source has been observed
- **url**: Link to AstroGeo database for detailed information

## Examples

### Example 1: Basic Phase Calibrator Search

```bash
planobs_phasecal -t 'M87' -b 6cm
```

This finds all phase calibrator candidates within 5 degrees of M87 at 6cm wavelength.

### Example 2: Tight Separation Constraint

```bash
planobs_phasecal -t 'M87' -b 6cm --max-separation 2 --min-flux 0.2
```

This finds calibrators within 2 degrees of M87 that have at least 0.2 Jy unresolved flux.

### Example 3: Multiple Band Search

```bash
planobs_phasecal -t '3C273' --max-separation 10 -n 10
```

This searches all bands for calibrators within 10 degrees of 3C273 and returns the top 10 closest candidates.

### Example 4: High Frequency Requirements

```bash
planobs_phasecal -t 'NGC 4258' -b 1.3cm --max-separation 3 --min-flux 0.5
```

This finds calibrators for a 1.3cm observation with stricter flux requirements suitable for high-frequency work.

## Phase Referencing Considerations

### Separation Distance

- **< 1 degree**: Excellent for phase calibration, minimal phase decorrelation
- **1-3 degrees**: Good, commonly used range
- **3-5 degrees**: Acceptable for lower frequencies
- **> 5 degrees**: May have significant phase decorrelation, especially at high frequencies

### Flux Requirements

- **Low frequencies (≥ 6cm)**: 0.1-0.2 Jy typically sufficient
- **Mid frequencies (3-6cm)**: 0.2-0.5 Jy recommended
- **High frequencies (≤ 2cm)**: 0.5-1.0+ Jy often needed

### Source Structure

- **Unresolved/Total flux ratio**: Higher ratios indicate more compact sources
- **Avoid extended sources**: Structure can introduce phase errors
- **Check AstroGeo links**: Look for information on source structure

## Tips for Good Phase Calibrators

1. **Proximity**: Choose the closest suitable calibrator
2. **Flux**: Ensure sufficient flux for your observing band and sensitivity
3. **Compactness**: Prefer sources with high unresolved/total flux ratios
4. **Elevation**: Consider calibrator elevation during your observation
5. **Multiple options**: Identify backup calibrators in case of issues

## Troubleshooting

### No Sources Found

If no sources are returned, try:

- Increasing the `--max-separation` distance
- Lowering the `--min-flux` threshold
- Checking if the target name is correct
- Trying a different observing band

### Target Not Found

If the target source is not recognized:

- Verify the target name spelling
- Try using coordinates instead of a name
- Check if the target is in the RFC catalog
- Use the main `planobs` tool to verify target coordinates

### Too Many Weak Sources

If you get many sources but they're all weak:

- Increase the `--min-flux` threshold
- Reduce the `--max-separation` distance
- Use the `-n` option to limit output to the best candidates

## Integration with Observations

Once you've identified suitable phase calibrators, you can:

1. **Add to source catalog**: Include them in your source catalog file
2. **Plan switching times**: Calculate slewing times between target and calibrators
3. **Design scan strategy**: Plan calibrator scans throughout your observation
4. **Verify with main tool**: Use the main `planobs` tool to verify visibility

The phase calibrator information integrates seamlessly with the main PlanObs workflow for complete observation planning.
