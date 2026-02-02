# Installation

## Requirements

- **Python 3.11+** (required)
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
planobs-server --help
```

## Available Commands

| Command | Description |
|---------|-------------|
| `planobs` | Command-line interface for observation planning |
| `planobs-server` | Launch the web-based GUI |

## Dependencies

PlanObs automatically installs these dependencies:

- **numpy** – Numerical computing
- **astropy** – Astronomical calculations
- **astroplan** – Observation planning
- **plotly/dash** – Interactive visualizations
- **rich** – Terminal formatting

## Troubleshooting

!!! warning "Python Version"
    PlanObs requires Python 3.11 or higher. Check your version with `python --version`.

!!! tip "Virtual Environment"
    We recommend using a virtual environment:
    ```bash
    python -m venv planobs-env
    source planobs-env/bin/activate  # Linux/macOS
    pip install vlbiplanobs
    ```

