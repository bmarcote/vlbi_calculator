# EVN Observation Planner (PlanObs)

A comprehensive tool for planning Very Long Baseline Interferometry (VLBI) observations.

---

## Overview

The **EVN Observation Planner** determines source visibility, estimates sensitivity, and calculates resolution for VLBI observations. While designed for the [European VLBI Network (EVN)](https://www.evlbi.org), it supports:

- [Very Long Baseline Array (VLBA)](https://public.nrao.edu/telescopes/vlba/)
- [Australian Long Baseline Array (LBA)](https://www.atnf.csiro.au/vlbi/overview/index.html)
- [eMERLIN](http://www.merlin.ac.uk/e-merlin/index.html)
- [Global mm-VLBI Array](https://www3.mpifr-bonn.mpg.de/div/vlbi/globalmm/)
- Custom ad-hoc VLBI arrays

!!! tip "Online Version"
    Use PlanObs directly at [planobs.jive.eu](https://planobs.jive.eu) — no installation required!

---

## Features

| Feature | Description |
|---------|-------------|
| **Visibility Plots** | See when sources are observable by each antenna |
| **Sensitivity Estimation** | Calculate expected RMS noise levels |
| **Resolution Calculation** | Estimate angular resolution for your array |
| **Scheduling** | Generate pySCHED-compatible `.key` files |
| **Multiple Interfaces** | GUI, CLI, and Python API |

---

## CLI Modes

The `planobs` command is organised into six modes:

| Mode | Command | Description |
|------|---------|-------------|
| **Observe** | `planobs -b 6cm -t 'M87' --network EVN` | Plan a VLBI observation |
| **Fringe Finders** | `planobs fringefinders -s Ef Hh Mc -e '2025-03-15 08:00' -d 8` | Find bright calibrators for fringe detection |
| **Phase Calibrators** | `planobs phasecals 'M87' -b 6cm` | Find compact calibrators near a target |
| **Source** | `planobs source '3C273'` | Look up source information |
| **Antenna** | `planobs antenna Ef` | Look up antenna information |
| **Server** | `planobs server` | Launch the web GUI |

=== "Python API"

    ```python
    from astropy.time import Time
    from astropy import units as u
    from vlbiplanobs import cli

    obs = cli.main(band='6cm', networks=['EVN'], targets=['J1230+1223'],
                   start_time=Time('2025-03-15 20:00', scale='utc'), duration=8*u.h)
    print(obs.thermal_noise())
    ```
    Full programmatic control for custom workflows.

---

## Quick Links

<div class="grid cards" markdown>

- :material-download: **[Installation](installation.md)** – Get PlanObs running locally
- :material-rocket-launch: **[Quick Start](quickstart.md)** – All modes in 5 minutes
- :material-telescope: **[Observation Planning](mode-observe.md)** – Full observe mode reference
- :material-satellite-antenna: **[Fringe Finders](fringefinder.md)** – Find bright calibrator sources
- :material-target-account: **[Phase Calibrators](phasecal.md)** – Find nearby phase calibrators
- :material-api: **[API Reference](api/index.md)** – Full Python API documentation

</div>

---

## About

PlanObs is developed at [JIVE](https://www.jive.eu) as the successor to the [EVN Calculator](http://old.evlbi.org/cgi-bin/EVNcalc.pl). It is open source under the GPL-3.0 license.

**Author:** Benito Marcote (marcote@jive.eu)

