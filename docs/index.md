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

## Usage Options

=== "GUI (Web Interface)"

    ```bash
    planobs-server
    ```
    Launches a graphical interface identical to the JIVE-hosted version.

=== "CLI (Command Line)"

    ```bash
    planobs -b 6cm -t 'M87' --network EVN
    ```
    Flexible command-line interface for scripting and automation.

=== "Python API"

    ```python
    from vlbiplanobs import Observation
    
    obs = Observation()
    obs.band = '6cm'
    obs.stations = obs.stations.filter_by_network('EVN')
    obs.add_target('M87')
    ```
    Full programmatic control for custom workflows.

---

## Quick Links

<div class="grid cards" markdown>

- :material-download: **[Installation](installation.md)** – Get PlanObs running locally
- :material-rocket-launch: **[Quick Start](quickstart.md)** – Your first observation in 5 minutes
- :material-book-open: **[Tutorial](tutorial.md)** – Detailed usage guide
- :material-api: **[API Reference](api/index.md)** – Full Python API documentation
- :material-satellite-antenna: **[Fringe Finders](fringefinder.md)** – Find bright calibrator sources
- :material-target-account: **[Phase Calibrators](phasecal.md)** – Find nearby phase calibrators

</div>

---

## About

PlanObs is developed at [JIVE](https://www.jive.eu) as the successor to the [EVN Calculator](http://old.evlbi.org/cgi-bin/EVNcalc.pl). It is open source under the GPL-3.0 license.

**Author:** Benito Marcote (marcote@jive.eu)

