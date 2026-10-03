# Installation

## Requirements

- **Python 3.12+** (required)
- pip or conda package manager

## Install from PyPI

The recommended way to install PlanObs:

```bash
pip install vlbiplanobs
```

## Install from Source

For development or the latest features:

```bash
git clone https://github.com/bmarcote/vlbi_calculator.git
cd vlbi_calculator
pip install -e .
```

## Verify Installation

After installation, verify the binaries are available:

```bash
planobs --help
planobs --version
planobs server --help
```

## Available Commands

After installation, the `planobs` command provides several modes:

| Command | Description |
|---------|-------------|
| `planobs` | Observation planning (default mode) |
| `planobs fringefinders` | Find fringe finder sources |
| `planobs phasecals` | Find phase calibrator sources |
| `planobs source` | Look up source information |
| `planobs antenna` (or `ant`) | Look up antenna information |
| `planobs server` | Launch the web-based GUI |

## Dependencies

PlanObs automatically installs these dependencies:

- **numpy** – Numerical computing
- **astropy**, **pyerfa** – Astronomical calculations (vectorized visibility computations)
- **astroplan** – Observation planning
- **ortools** – Optimization in the scheduler
- **plotly/dash** (with dash-bootstrap-components and dash-mantine-components) – Web GUI and interactive plots
- **kaleido** (< 1) and **borb** – PDF reports (kaleido ≥ 1 needs a Chrome binary, so it is pinned below 1)
- **rich**, **rich_argparse**, **plotext** – Terminal output and plots

## Troubleshooting

!!! warning "Python Version"
    PlanObs requires Python 3.12 or higher. Check your version with `python --version`.

!!! tip "Virtual Environment"
    We recommend using a virtual environment:
    ```bash
    python -m venv planobs-env
    source planobs-env/bin/activate  # Linux/macOS
    pip install vlbiplanobs
    ```

