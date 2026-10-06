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
| `--logging [LOGFILE]` | off | Log to a file (default `/var/log/planobs.log` if writable, otherwise `~/log-planobs.log`). |

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

## Environment variables

| Variable | Default | Description |
|----------|---------|-------------|
| `PLANOBS_POLARIS_ORIGINS` | unset | Comma-separated list of web origins (e.g. `https://polaris.example.org,https://test.example.org:8443`) allowed to receive the observation setup via **Export to Polaris**. When set, the export is offered only if PlanObs was opened from one of these origins; otherwise the button falls back to **Copy link**. When unset, the opener origin taken from the page referrer is used. The setup is always posted to that exact origin, never to `*`. |

---

## Production deployment

`planobs server` runs the Dash development server. For a public deployment, serve the WSGI app
`vlbiplanobs.gui.main:server` with a production server such as gunicorn, using the settings shipped with the
package (this is the command used in the `Procfile` of the repository):

```bash
gunicorn -c python:vlbiplanobs.gui.gunicorn_conf vlbiplanobs.gui.main:server
```

These settings preload the app in the gunicorn master process (`preload_app`), so that the warm-up done when the
app is imported (loading the Earth-rotation table, the astropy/astroplan/plotly internals, and the page layout)
happens only once and every worker answers its very first request at full speed. They can be tuned with
environment variables:

| Variable | Default | Description |
|----------|---------|-------------|
| `PLANOBS_BIND` | `127.0.0.1:8050` | Address and port to listen on. |
| `PLANOBS_WORKERS` | `4` (or the number of CPUs if lower) | Number of worker processes. |
| `PLANOBS_THREADS` | `4` | Threads per worker. |
| `PLANOBS_TIMEOUT` | `120` | Seconds before a busy worker is restarted (PDF reports can take several seconds). |
| `PLANOBS_MAX_REQUESTS` | `5000` | Requests served by a worker before it is replaced. |
| `PLANOBS_NO_WARMUP` | unset | Set it to skip the warm-up at import time (faster startup, slower first request). |

A systemd unit would then contain:

```ini
[Service]
User=www-data
Group=www-data
WorkingDirectory=/var/local/gunicorn/
ExecStart=/path/to/env/bin/gunicorn -c python:vlbiplanobs.gui.gunicorn_conf vlbiplanobs.gui.main:server
```

The app compresses its responses when `flask-compress` is installed (it is a dependency, through
`dash[compress]`) and lets the browsers cache the static files. When running behind nginx, it is more
efficient to let nginx do both, e.g.:

```nginx
gzip on;
gzip_types application/json application/javascript text/css image/svg+xml;
gzip_min_length 1024;

location /assets/ {
    alias /path/to/env/lib/python3.x/site-packages/vlbiplanobs/gui/assets/;
    expires 1d;
}
location / {
    proxy_pass http://127.0.0.1:8050/;
    include proxy_params;
}
```

Station and network catalogs are shared by all users of a worker process, but each user's observation
(including per-station data rates) is computed independently.

---

## Input limits and validation

All inputs are validated on the server, including those coming from shared `?config=` links, the browser's
local storage, and uploaded files. Invalid values are rejected with an explanatory message. Current limits:

| Input | Limit |
|-------|-------|
| Observation duration | ≤ 50 h |
| Number of target sources | ≤ 20 |
| Source name length | ≤ 80 characters |
| Uploaded source-list file | ≤ 50 kB |
| Setup values (band, data rate, subbands, channels, polarizations, integration time) | Only the values offered in the GUI |

Source names sent to the online resolvers (SIMBAD/NED/VizieR) are validated first and the results are cached.
In the generated SCHED `.key` files, quotes and newlines are removed from all user-provided strings.

No analytics or tracking is used.

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
