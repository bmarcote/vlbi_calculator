# Web Server (GUI)

The **server** mode launches the PlanObs web interface — the same application hosted at [planobs.jive.eu](https://planobs.jive.eu). It runs a local Dash/Plotly server you can access in your browser.

## Usage

```bash
planobs server [OPTIONS]
```

---

## Options

| Argument | Default | Description |
|----------|---------|-------------|
| `--host` | `127.0.0.1` | Network interface to bind to. Use `0.0.0.0` to allow access from other machines. |
| `--port` | `8050` | TCP port number. |
| `--debug` | off | Enable Dash debug mode (auto-reload on code changes, detailed error pages). |

---

## Examples

### Start with defaults

```bash
planobs server
```

Then open [http://localhost:8050](http://localhost:8050) in your browser.

### Expose on the local network

```bash
planobs server --host 0.0.0.0 --port 8080
```

Other machines on the same network can access the GUI at `http://<your-ip>:8080`.

### Development mode

```bash
planobs server --debug
```

The server auto-reloads when source files change and shows detailed tracebacks on errors.

---

## Features

The web GUI provides interactive versions of the same capabilities available through the CLI:

- **Source visibility plots** – interactive elevation plots for each station.
- **Sensitivity calculator** – expected thermal noise based on the array and setup.
- **Resolution estimator** – angular resolution for the selected stations and band.
- **Station map** – geographic map of the participating antennas.

### Antenna Selection

The antenna panel lets you select stations individually or by network.

- **Single-configuration antennas** appear as standard toggle chips. Click to include/exclude from the observation.
- **Multi-configuration antennas** (VLA, MeerKAT, SKAO) appear as a **split-button chip**:
  - The chip label shows the group name (e.g. `VLA`).
  - The **▼** button on the right opens a dropdown to switch between configurations (e.g. single dish vs. phased array). Hovering over each option shows the full antenna information card.
  - Selecting any configuration from the dropdown automatically includes that antenna in the observation.
  - The currently active configuration is highlighted in the dropdown menu.
  - Only one configuration per group can be active at a time.
- Chips dim automatically when the active configuration does not support the selected observing band.
- **Network switches** select/deselect all antennas in a predefined network at once, including the appropriate configuration for grouped antennas.

---

## When to Use This Mode

- **Interactive exploration** – experiment with different arrays, bands, and times without re-running CLI commands.
- **Presentations and teaching** – the graphical output is suitable for talks and tutorials.
- **Remote access** – bind to `0.0.0.0` to share the planner with colleagues on your network.
