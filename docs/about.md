# About PlanObs

## Overview

The EVN Observation Planner (PlanObs) is a comprehensive tool for planning Very Long Baseline Interferometry (VLBI) observations. It was developed at [JIVE](https://www.jive.eu) as the successor to the EVN Calculator.

## Features

- **Source Visibility** – Determine when astronomical sources are visible by different antennas
- **Sensitivity Estimation** – Calculate expected RMS noise levels for planned observations
- **Resolution Calculation** – Estimate the angular resolution of your VLBI array
- **Schedule Generation** – Create pySCHED-compatible `.key` files for EVN observations
- **Multiple Interfaces** – GUI, CLI, and Python API for different workflows

## Supported Networks

PlanObs supports all major VLBI networks:

- European VLBI Network (EVN)
- Very Long Baseline Array (VLBA)
- Australian Long Baseline Array (LBA)
- eMERLIN
- Korean VLBI Network (KVN)
- Global mm-VLBI Array (GMVA)
- Custom ad-hoc arrays

## Development

PlanObs is developed and maintained by Benito Marcote at JIVE.

- **Repository:** [github.com/bmarcote/vlbi_calculator](https://github.com/bmarcote/vlbi_calculator)
- **Issues:** [GitHub Issues](https://github.com/bmarcote/vlbi_calculator/issues)
- **Contact:** marcote@jive.eu

## Acknowledgments

PlanObs uses several open-source libraries:

- [Astropy](https://www.astropy.org/) – Astronomical calculations
- [Astroplan](https://astroplan.readthedocs.io/) – Observation planning
- [Plotly/Dash](https://dash.plotly.com/) – Interactive visualizations
- [NumPy](https://numpy.org/) – Numerical computing

## Citation

If you use PlanObs in your research, please cite:

```
Marcote, B. (2024). EVN Observation Planner (PlanObs).
https://github.com/bmarcote/vlbi_calculator
```
