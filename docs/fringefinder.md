# Fringe Finder Sources

The `planobs_fringefinder` tool helps you find suitable fringe finder sources for VLBI observations. Fringe finders are bright calibrator sources used to detect and correct instrumental delays between telescopes.

## Overview

Fringe finder sources must be:
- Bright enough to produce strong fringes across all baselines
- Visible by all participating stations during the observation
- Above minimum elevation requirements
- Located at suitable positions for calibration

## Usage

### Basic Usage

```bash
planobs_fringefinder -s Ef Hh Mc Tr -t '2025-03-15 08:00' -d 8 -b 6cm
```

### Advanced Usage

```bash
planobs_fringefinder -s Ef Hh Mc Tr Ys -t '2025-03-15 08:00' -d 8 -b 6cm \
    --min-flux 1.0 --min-elevation 25 --require-all -l 10
```

## Command Line Options

### Required Options

- `-s, --stations`: List of antenna codenames that will participate in the observation
  - Example: `Ef Hh Mc Tr` or `Effelsberg Hartebeesthoek Medicina Noto`
  - Use station codenames as defined in the station catalog

- `-t, --starttime`: Start of the observation
  - Format: 'YYYY-MM-DD HH:MM' (UTC)
  - Example: '2025-03-15 08:00'

- `-d, --duration`: Duration of the observation in hours
  - Example: `8` for an 8-hour observation

- `-b, --band`: Observing band
  - Supported bands: '18cm', '21cm', '13cm', '6cm', '5cm', '3.6cm', '2cm', '1.3cm', '0.7cm'
  - The tool will automatically convert wavelength bands to RFC catalog bands

### Optional Options

- `--min-flux`: Minimum unresolved flux threshold in Jy (default: 0.5)
  - Higher values return fewer, brighter sources
  - Lower values return more candidates but may include weaker sources

- `--min-elevation`: Minimum elevation in degrees (default: 20)
  - Sources must be above this elevation for all stations
  - Higher values reduce the number of visible sources

- `-l, --max-lines`: Maximum number of sources to return (default: 20)
  - Limits the output to the top N candidates
  - Sources are sorted by flux density

- `--require-all`: Require source to be visible by ALL stations (default: False)
  - When set, sources must be visible by every specified station
  - When not set, sources visible by any station are included

- `--station-catalog`: Path to custom station catalog file
  - Use if you have a modified station catalog
  - Overrides the default catalog

## Output

The tool displays a table with the following columns:

- **Name**: Source name from the RFC catalog
- **IVS Name**: International VLBI Service name
- **Min elev. (deg)**: Minimum elevation across all stations
- **Total flux (Jy)**: Total flux density at the observing band
- **Unresolved (Jy)**: Unresolved flux density (most relevant for fringe finding)
- **Bands**: Bands where the source has been observed
- **url**: Link to AstroGeo database for detailed information

## Examples

### Example 1: EVN Fringe Finders

```bash
planobs_fringefinder -s Ef Hh Mc Tr Wb On -t '2025-06-15 20:00' -d 12 -b 6cm
```

This finds fringe finder candidates for a 12-hour EVN observation at 6cm, requiring visibility from Effelsberg, Hartebeesthoek, Medicina, Noto, Wettzell, and Onsala.

### Example 2: High Elevation Requirement

```bash
planobs_fringefinder -s Ef Hh Mc Tr -t '2025-06-15 20:00' -d 8 -b 6cm \
    --min-elevation 30 --min-flux 0.8
```

This finds sources that are above 30 degrees elevation and have at least 0.8 Jy unresolved flux.

### Example 3: Strict Requirements

```bash
planobs_fringefinder -s Ef Hh Mc Tr Ys -t '2025-06-15 20:00' -d 8 -b 6cm \
    --require-all --min-flux 1.0 -l 5
```

This finds only the top 5 sources that are visible by ALL stations and have >1.0 Jy flux.

## Tips for Good Fringe Finders

1. **Flux**: Look for sources with high unresolved flux (>0.5 Jy typically)
2. **Elevation**: Higher elevation sources are preferred (less atmospheric effects)
3. **Visibility**: Sources should be visible for the entire observation period
4. **Position**: Sources distributed across the sky can help with calibration
5. **Structure**: Point-like sources (high unresolved/total flux ratio) are ideal

## Troubleshooting

### No Sources Found

If no sources are returned, try:
- Lowering the `--min-flux` threshold
- Reducing the `--min-elevation` requirement
- Removing the `--require-all` flag
- Checking if all stations can observe at the specified band

### Too Many Sources

If too many sources are returned, try:
- Increasing the `--min-flux` threshold
- Using the `-l` option to limit output
- Adding the `--require-all` flag for stricter requirements

### Station Not Found

If a station name is not recognized:
- Check the station codename in the station catalog
- Use the `--list-antennas` option with the main `planobs` command to see available stations
