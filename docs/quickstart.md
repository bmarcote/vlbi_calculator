# Quick Start

Get started with PlanObs in under 5 minutes. The `planobs` command provides five modes, each designed for a specific task in the VLBI observation planning workflow.

## Modes Overview

| Mode | Command | Purpose |
|------|---------|---------|
| **Observe** | `planobs -b BAND ...` or `planobs observe ...` | Plan a VLBI observation: visibility, sensitivity, scheduling |
| **Fringe Finders** | `planobs fringefinders ...` | Find bright calibrator sources for fringe detection |
| **Phase Calibrators** | `planobs phasecals ...` | Find compact calibrators near your target |
| **Source** | `planobs source <name>` | Look up detailed information about a source |
| **Server** | `planobs server` | Launch the web-based GUI |

!!! tip "Legacy syntax"
    Running `planobs` without a subcommand (e.g. `planobs -b 6cm ...`) is equivalent to `planobs observe -b 6cm ...`. Both forms work.

---

## 1. Observe – Plan a VLBI Observation

Check when a source is visible with the EVN and estimate the expected sensitivity:

```bash
planobs -b 6cm -t 'M87' --network EVN
```

This prints source visibility per antenna, expected thermal noise, and optimal observing windows.

:material-arrow-right: **[Full observe reference](mode-observe.md)**

---

## 2. Fringe Finders – Find Calibrator Sources

Find bright, compact sources visible by your stations during the observation, suitable for fringe detection:

```bash
planobs fringefinders -s Ef Hh Mc Tr -t '2025-03-15 08:00' -d 8 -b 6cm
```

Returns a ranked table of candidates with flux densities, elevations, and AstroGeo links.

:material-arrow-right: **[Full fringefinders reference](fringefinder.md)**

---

## 3. Phase Calibrators – Find Nearby Calibrators

Search for phase calibrator candidates close to your target:

```bash
planobs phasecals -t 'M87' -b 6cm
```

Returns calibrators within the default 5° separation, sorted by distance.

:material-arrow-right: **[Full phasecals reference](phasecal.md)**

---

## 4. Source – Look Up Source Information

Retrieve RFC catalog data (flux densities, coordinates, AstroGeo link) for any known source:

```bash
planobs source '3C273'
```

If the source is not in the RFC catalog, PlanObs resolves it via SIMBAD/NED/VizieR.

:material-arrow-right: **[Full source reference](mode-source.md)**

---

## 5. Server – Launch the Web GUI

Start the interactive web interface (identical to [planobs.jive.eu](https://planobs.jive.eu)):

```bash
planobs server
```

Then open your browser to `http://localhost:8050`.

:material-arrow-right: **[Full server reference](mode-server.md)**

---

## Next Steps

- **[Observation Planning](mode-observe.md)** – Complete guide to the observe mode
- **[Scheduling](scheduling.md)** – Generate pySCHED `.key` files
- **[CLI Overview](cli.md)** – All modes and common patterns at a glance
- **[API Reference](api/index.md)** – Python API documentation
