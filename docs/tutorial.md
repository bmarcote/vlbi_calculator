# Tutorial

A step-by-step walkthrough covering the most common PlanObs workflows.

!!! tip
    For a quick overview of every mode with one-liner examples, see the **[Quick Start](quickstart.md)**.

---

## 1. Plan a VLBI Observation

### Check source visibility (no fixed time)

Find when a source is observable by the EVN and eMERLIN:

```bash
planobs -b 6cm -t 'Altair' --network EVN eMERLIN
```

PlanObs scans the full GST range and reports optimal observing windows.

### Plan at a specific epoch

```bash
planobs -b 6cm -t 'Altair' --stations Ef Hh Ir Mc Tr T6 O8 Wb Cm \
  --starttime '2025-06-15 20:00' --duration 12
```

This evaluates visibility and sensitivity for the given 12-hour window.

### Suppress the GUI

Add `--no-tui` to skip the graphical elevation plots and keep only the terminal output:

```bash
planobs -b 6cm -t 'Altair' --network EVN --no-tui
```

:material-arrow-right: Full reference: **[Observation Planning](mode-observe.md)**

---

## 2. Find Fringe Finder Sources

Fringe finders are bright calibrators used to detect and correct instrumental delays. Use the `fringefinders` mode:

```bash
planobs fringefinders -s Ef Hh Mc Tr -t '2025-03-15 08:00' -d 8 -b 6cm
```

This searches the RFC catalog and returns a ranked table of candidates with flux densities, elevations, and AstroGeo links.

:material-arrow-right: Full reference: **[Fringe Finders](fringefinder.md)**

---

## 3. Find Phase Calibrators

Phase calibrators are compact sources near the target used for phase referencing:

```bash
planobs phasecals -t 'M87' -b 6cm --max-separation 5 --min-flux 0.1
```

Returns calibrators within 5° of M87 with ≥ 0.1 Jy unresolved flux.

:material-arrow-right: Full reference: **[Phase Calibrators](phasecal.md)**

---

## 4. Look Up Source Information

Query the RFC catalog (or SIMBAD/NED/VizieR) for any source:

```bash
planobs source '3C273'
```

Prints coordinates, flux table across bands, number of observations, and the AstroGeo link.

:material-arrow-right: Full reference: **[Source Lookup](mode-source.md)**

---

## 5. Generate a Schedule File

Create a pySCHED-compatible `.key` file by adding `--sched`:

```bash
planobs -b 6cm -t 'M87' --network EVN \
  --starttime '2025-03-15 08:00' --duration 8 \
  --sched eg123a
```

This produces `eg123a.key`. You can control fringe finder selection and add calibration scans:

```bash
planobs -b 6cm -t 'M87' --network EVN \
  --starttime '2025-03-15 08:00' --duration 8 \
  --sched eg123a --fringefinders 3 --polcal
```

:material-arrow-right: Full reference: **[Scheduling](scheduling.md)**

---

## 6. Customising the Station Catalog

Station data is stored in `data/stations_catalog.inp` using [Python configparser format](https://docs.python.org/3/library/configparser.html). To add a station, create a new section:

```ini
[Station Name]
station = Station Name
code = XX
network = EVN
possible_networks = EVN, GMVA
country = Netherlands
diameter = 25 m
position = X, Y, Z
min_elevation = 10
real_time = yes
SEFD_6 = 420
SEFD_18 = 700
```

| Field | Description |
|-------|-------------|
| `code` | Unique codename (e.g. `Ef`, `Wb`). |
| `network` | Primary network, or `Other`. |
| `possible_networks` | Comma-separated list of networks the station can join. |
| `position` | Geocentric X, Y, Z coordinates. |
| `min_elevation` | Minimum observable elevation in degrees (default 10). |
| `real_time` | `yes` if the station supports e-VLBI. |
| `SEFD_YY` | System Equivalent Flux Density at wavelength YY cm (in Jy). One entry per band. |

Use `--station-catalog` to load your custom file:

```bash
planobs -b 6cm --station-catalog my_stations.inp --network EVN -t 'M87'
```

We welcome contributions of additional stations — please open an issue or pull request on [GitHub](https://github.com/bmarcote/vlbi_calculator).

