# Source Lookup

The **source** mode retrieves detailed information about a single source. It queries the RFC (Radio Fundamental Catalog) and, if the source is not found there, falls back to SIMBAD/NED/VizieR.

## Usage

```bash
planobs source [options] <source_name>
```

The `source_name` argument is **required** and positional. It can be:

- A J2000 name from the RFC catalog (e.g. `J1229+0203`)
- An IVS name (e.g. `3C273`)
- Any source resolvable by SIMBAD/NED/VizieR (e.g. `M87`)

### Options

| Flag | Description |
|------|-------------|
| `--no-networks` | Skip the network observability table. |
| `--gst` | Append a **GST** column to the observability table showing the Greenwich Sidereal Time ranges (`HH:MM-HH:MM`) during which more than 3 antennas in the network can observe the source. |

---

## Output

### Source found in RFC catalog

When the source is in the RFC catalog, PlanObs prints:

- **Name and aliases** – J2000 name and IVS name.
- **Coordinates** – Right ascension and declination in HMS/DMS format.
- **Number of observations** – how many times the source has been observed in geodetic/astrometric VLBI.
- **Flux table** – total and unresolved flux density at each band where data is available.
- **AstroGeo link** – direct URL to the source page on the AstroGeo database.

The flux table includes these bands:

| Band | Wavelength |
|------|------------|
| S | 18/21 cm |
| C | 13/6/5 cm |
| X | 3.6 cm |
| U | 2 cm |
| K | 1.3/0.7 cm |

### Source not in RFC catalog

If the source is not found in the RFC catalog but can be resolved by SIMBAD/NED/VizieR, PlanObs prints:

- **Name** – as returned by the resolver.
- **Coordinates** – right ascension and declination.
- A note that the source is not in the RFC calibrators catalog.

If the source cannot be found at all, PlanObs prints an error and exits with code 1.

---

## Examples

### Look up a well-known calibrator

```bash
planobs source '3C273'
```

### Look up by J2000 name

```bash
planobs source 'J1229+0203'
```

### Look up a non-RFC source

```bash
planobs source 'M87'
```

### Show GST observing windows

```bash
planobs source '3C273' --gst
```

The table gains a **GST** column showing the sidereal time range(s) when each network has more than 3 antennas able to observe the source (e.g. `01:15-19:30`).

---

## When to Use This Mode

- **Before planning an observation** – check whether a potential calibrator has sufficient flux.
- **Evaluating calibrators** – compare total vs. unresolved flux to judge compactness.
- **Quick coordinate lookup** – get precise J2000 coordinates for any source.
