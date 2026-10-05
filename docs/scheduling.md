# Scheduling

PlanObs includes a scheduler that generates pySCHED-compatible `.key` files for EVN observations.

## Overview

The scheduler arranges scan blocks across your observation following VLBI conventions:

- **Fringe finders**: Auto-selected or specified sources, distributed throughout the observation
- **Phase calibrators**: Auto-selected or specified sources for phase referencing
- **Check sources**: Auto-selected or specified sources for calibration verification
- **Polarization calibration**: Standard polcal sources (3C84, OQ208, DA193) near 10%, 50%, 90% of observation time,
  at the closest time when each source is visible (polcals that cannot be placed are skipped with a warning)
- **eMERLIN 3C286**: Automatically added when any of the eMERLIN stations Cm, Da, De, Kn or Pi is present (Jodrell
  Bank, Jb1/Jb2, does not trigger it as it observes regularly in the EVN), at a time when it is visible
- **Science blocks**: Optimized for antenna participation and elevation

## Generating a Schedule File

### Via CLI

```bash
planobs -b 6cm -t 'M87' --network EVN \
  --epoch '2025-03-15 08:00' \
  --duration 8 \
  --sched my_experiment
```

This creates `my_experiment.key`.

### Frequency setup

Use `--setup` to write a given frequency setup in the `setup = ...` line of the `.key` file:

```bash
planobs -b 6cm -t 'M87' --network EVN \
  --epoch '2025-03-15 08:00' \
  --duration 8 \
  --sched my_experiment \
  --setup 'evn6cm-2Gbps-32MHz.set'
```

If `--setup` is not given, PlanObs guesses the setup from the observation, and writes a
`nosetup` placeholder in the `.key` file when it cannot determine one.

### Via Python

```python
from astropy.time import Time
from astropy import units as u
from vlbiplanobs import cli
from vlbiplanobs.scheduler import ObservationScheduler

obs = cli.main(band='6cm', networks=['EVN'], targets=['J1230+1223'],
               start_time=Time('2025-03-15 20:00', scale='utc'), duration=8*u.h)

scheduler = ObservationScheduler(obs, fringefinder_spec=['3'], polcal=True)
scheduler.schedule()   # must be called before generate_key_file()
key_content = scheduler.generate_key_file(experiment_code='EG123A', pi_name='Your Name',
                                          pi_email='you@example.com', pi_institute='Your Institute',
                                          setup_file='evn6cm-2Gbps-32MHz.set')

with open('eg123a.key', 'w') as f:
    f.write(key_content)
```

User-provided strings (experiment code, PI fields, source names) have quotes and newlines removed before they
are written. With `--sched`, the experiment code is the file name stem in upper case and may only contain
letters, digits and underscores.

## Schedule File Format

The generated `.key` file follows the pySCHED format:

```
! =================  Cover Information  ====================
expcode  = 'EG123A'
piname   = 'PI Name'
email    = 'pi@example.com'
obstype  = 'VLBI'

! ==============  Correlator Information  ==================
correl   = 'JIVE'
cornant  = '13'

srccat /
  source='J1224+2122' ra=12:24:54.4584 dec=+21:22:46.388 equinox='J2000' /
  source='J1230+1223' ra=12:30:49.4234 dec=+12:23:28.044 equinox='J2000' /
endcat /

! ==================  Frequency Setup  =====================
setup = 'evn6cm-2Gbps-32MHz.set'

! ===============  Start of Observation  ===================
year  = 2025
month = 3
day   = 15
start = 08:00:00

! =====================  The Scans  ========================
stations = eflsberg, hart, jodrell2, medicina, noto, onsala85, torun, yebes40m, wstrbork

source='J1224+2122' gap=0:00 dur=5:00 intent='FRINGE_FINDER' /
source='J1225+1253' gap=0:00 dur=1:30 intent='PHASE_CAL' /
source='J1230+1223' gap=0:00 dur=3:30 intent='TARGET' /
...
```

This is an abridged example. Stations are written with their pySCHED names, and with `--sched` the PI fields keep
the template placeholders (`PI Name`, `pi@example.com`), to be edited by hand.

### Custom template

`--template FILE` uses your own `.key` template instead of the bundled one. Get a copy of the bundled template to
start from with:

```bash
planobs --get-key-template my_template.key
```

## Using with pySCHED

After generating the `.key` file:

```bash
sched.py -k my_experiment.key
```

This produces the VEX file and other outputs needed for correlation.

## Scheduler Configuration

### Auto-Selection Features

The scheduler can automatically select calibrators for your observation:

- **Fringe finders**: Use `--fringefinders N` to schedule N fringe-finder scans (default 2) on a bright calibrator auto-selected from the RFC catalog, or specify named sources with `--fringefinders '3C273' '3C279'`.
- **Phase calibrators**: Use `--phasecal` (empty) to auto-select the best phase calibrator based on unresolved flux, compactness, and proximity to the target. Specify named sources with `--phasecal 'J1229+0203'`.
- **Check sources**: Use `--check-source` (empty) to auto-select a check source close to the target with similar properties to the phase calibrator.

### Fringe Finder Rules

| `--fringefinders` | Fringe finder scans |
|-------------------|---------------------|
| not given | 2 scans on one auto-selected source |
| a number `N` | `N` scans on one auto-selected source |
| source names | 1 scan up to 1.5 h, 2 up to 3 h, and one more per additional ~2 h |

Each fringe finder scan is 5 minutes, and the scans are spread across the observation. When multiple fringe finders are specified, they are distributed round-robin throughout the observation.

### Science Block Optimization

The scheduler optimizes science blocks for:

1. **Maximum antennas** – Prioritizes times when most stations can observe
2. **Highest elevation** – Secondary optimization for better sensitivity

### Python API

```python
scheduler = ObservationScheduler(
    obs,
    min_antennas=3,           # Require at least 3 antennas
    require_all_antennas=True # All stations must observe
)

schedule = scheduler.schedule()
for block in scheduler.get_scheduled_blocks():
    print(block.name, block.start_time.iso, block.end_time.iso, block.n_antennas)
```

## Network Monitoring Experiments

`planobs ... --nme [--sched CODE]` produces NME schedules (all-antenna fringe-finder scans with periodic ftp
fringe-test grabs). See **[Observation Planning](mode-observe.md#network-monitoring-experiments-nme)** for the rules,
and `vlbiplanobs.nme` for the Python API (`plan_nme`, `generate_nme_key_file`).

## References

- [EVN Schedule Preparation](https://www.evlbi.org/observing-schedule-preparation)
- [NRAO SCHED Manual](https://science.nrao.edu/facilities/vlba/docs/manuals/obsvlba/schedfile)
