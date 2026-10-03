# -*- coding: utf-8 -*-
# Licensed under GPLv3+ - see LICENSE
"""VLBI Observation Scheduler.

Arranges scan blocks across a VLBI observation: fringe finders, polarisation
calibrators, eMERLIN 3C286, and science targets. Generates SCHED .key files
with ``group N rep R`` syntax and Jb1 source-change mitigation.

Classes
-------
ScheduledScanBlock
    Dataclass for a scan block placed at a specific time.
ObservationScheduler
    Main scheduler that orchestrates placement and key-file generation.
"""

from dataclasses import dataclass, field
from typing import Optional
from importlib import resources
from datetime import datetime
import logging
import os
import re
from pathlib import Path
import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from .sources import Source, SourceType, ScanBlock, Scan
from .observation import Observation
from . import calibrators

try:
    from ortools.sat.python import cp_model
    _HAS_ORTOOLS = True
except ImportError:  # pragma: no cover - exercised only when ortools missing
    cp_model = None  # type: ignore
    _HAS_ORTOOLS = False

log = logging.getLogger(__name__)

# eMERLIN stations that trigger the 3C286 flux-scale scan. Jb2 is excluded on purpose: it observes regularly
# within the EVN, so an EVN-only array with Jb2 must not get an eMERLIN 3C286 scan.
_EMERLIN_CODES = {'CM', 'KN', 'PI', 'DA', 'DE', 'JB1'}
_POLCAL_NAMES = ['3C84', 'OQ208', 'DA193']
_3C286_COORD = '13h31m08.288s +30d30m32.96s'
_JB1_MAX_SRC_CHANGES_PER_HOUR = 12


def _fmt_dur(q: u.Quantity) -> str:
    """Format an astropy duration quantity as 'M:SS' for SCHED key files.

    Parameters
    ----------
    q : Quantity
        Duration quantity.

    Returns
    -------
    str
        Formatted duration string.
    """
    total_sec = int(round(q.to(u.s).value))
    return f"{total_sec // 60}:{total_sec % 60:02d}"


def _sched_safe(text: object) -> str:
    """Make a user-provided string safe to insert into a SCHED single-quoted field or cover letter.

    SCHED values are single-quoted and the key file is line-oriented, so a newline or a single quote in
    user input could terminate the field and inject arbitrary SCHED commands.

    Parameters
    ----------
    text : object
        Value to sanitise (converted with ``str``).

    Returns
    -------
    str
        The text with CR/LF replaced by spaces and single quotes removed.
    """
    return str(text).replace('\r\n', ' ').replace('\r', ' ').replace('\n', ' ').replace("'", '')


def format_setup_line(setup_file: str | None) -> str:
    """Build the value for the frequency-setup line of a SCHED key file.

    The returned string is meant to be inserted into a template line such as
    ``setup = {SETUP}``.

    Parameters
    ----------
    setup_file : str or None
        Frequency setup given by the user (e.g. 'evn6cm-2Gbps-32MHz.set' or 'EFF_BAND_32').
        Surrounding quotes and whitespace are ignored. If None or empty, a
        'nosetup' placeholder is returned instead.

    Returns
    -------
    str
        Either "'<setup_file>'" or "nosetup   ! TODO: Add frequency setup".
    """
    if setup_file is None:
        return "nosetup   ! TODO: Add frequency setup"

    cleaned = setup_file.strip().strip('\'"').strip()
    if not cleaned:
        return "nosetup   ! TODO: Add frequency setup"

    return f"'{cleaned}'"


def _pysched_setup_dir() -> Path:
    """Return the pySCHED setups directory, creating the path if absent."""
    return Path.home() / '.pysched' / 'setups'


def _canonical_band(band: str) -> str:
    """Return a canonical band string used to match setup file names.

    Strips trailing whitespace and lowercases the input, preserving digits
    and the 'cm'/'mm' suffix (e.g. '18cm' -> '18cm', '5cm' -> '5cm').
    """
    return band.strip().lower()


def _parse_setup_filename(filename: str) -> dict[str, object]:
    """Parse a pySCHED setup filename into searchable features.

    Recognises patterns such as ``evn18cm-2Gbps-32MHz.set``,
    ``lba5cm-2p-2IF.set`` and VLBA-style ``v18cm-512-16-2.set``.

    Parameters
    ----------
    filename : str
        Name of the setup file (with or without path).

    Returns
    -------
    dict
        Dictionary with ``band``, ``rate_mbps``, ``pols``, ``ifs``,
        ``chan_bw_mhz`` and ``prefix`` keys.  Unmatched fields are None.
    """
    name = os.path.basename(filename).lower()
    base, _ = os.path.splitext(name)
    result: dict[str, object] = {
        'filename': filename, 'base': base, 'band': None, 'rate_mbps': None,
        'pols': None, 'ifs': None, 'chan_bw_mhz': None, 'prefix': None,
    }

    # Band: digits optionally with a decimal point, followed by cm or mm.
    band_match = re.search(r'(\d+(?:\.\d+)?)(cm|mm)', base)
    if band_match:
        result['band'] = band_match.group(1) + band_match.group(2)
        result['prefix'] = base[:band_match.start()].lower()

    # Data rate: 2Gbps, 512Mbps, 1G, etc.
    rate_match = re.search(r'(\d+(?:\.\d+)?)\s*(gbps|g|mbps|m)', base)
    if rate_match:
        value = float(rate_match.group(1))
        unit = rate_match.group(2).lower()
        result['rate_mbps'] = value * 1e3 if unit.startswith('g') else value

    # LBA-style polarizations and IFs.
    pol_match = re.search(r'(\d+)p', base)
    if pol_match:
        result['pols'] = int(pol_match.group(1))
    if_match = re.search(r'(\d+)if', base)
    if if_match:
        result['ifs'] = int(if_match.group(1))

    # EVN-style channel bandwidth: e.g. 32MHz, 16MHz.
    cbw_match = re.search(r'(\d+(?:\.\d+)?)\s*mhz', base)
    if cbw_match:
        result['chan_bw_mhz'] = float(cbw_match.group(1))

    return result


def _network_prefix(observation: Observation) -> str:
    """Infer the preferred setup-file prefix from the station composition.

    Uses ``Observation.guess_network`` to rank networks by antenna overlap and
    returns a conventional setup-file prefix ('evn', 'lba', 'v' for VLBA, or
    lowercased network name for others).
    """
    networks = Observation.guess_network(observation.band, list(observation.stations))
    if not networks:
        return ''
    primary = networks[0].lower()
    mapping = {'evn': 'evn', 'lba': 'lba', 'vlba': 'v'}
    return mapping.get(primary, primary)


def _is_global(observation: Observation) -> bool:
    """Return True if non-EVN VLBA stations dominate the array."""
    networks = Observation.guess_network(observation.band, list(observation.stations))
    if not networks:
        return False
    return networks[0].upper() != 'EVN'


def _setup_file_score(parsed: dict[str, object], observation: Observation,
                      preferred_prefix: str) -> float:
    """Score a candidate setup file against the observation parameters.

    Higher scores indicate better matches.  A negative score means the file is
    unsuitable (wrong band or array prefix).
    """
    if parsed['band'] != _canonical_band(observation.band):
        return -1.0

    prefix = parsed.get('prefix') or ''
    prefix_lower = prefix.lower()

    # Prefix must match the array family.
    if preferred_prefix == 'evn' and not prefix_lower.startswith('evn'):
        return -1.0
    if preferred_prefix == 'lba' and not prefix_lower.startswith('lba'):
        return -1.0
    if preferred_prefix == 'v' and not (prefix_lower.startswith('v') or prefix_lower.startswith('vlba')):
        return -1.0
    if preferred_prefix and not prefix_lower.startswith(preferred_prefix):
        return -1.0

    score = 0.0

    # Data rate is the strongest discriminator.
    datarate_mbps = observation.datarate.to(u.Mbit / u.s).value if observation.datarate is not None else None
    if datarate_mbps is not None and parsed['rate_mbps'] is not None:
        score += 100.0 - abs(datarate_mbps - parsed['rate_mbps']) / max(datarate_mbps, 1.0)

    # Prefer global EVN setup files when VLBA stations are present.
    if preferred_prefix == 'evn' and _is_global(observation) and '+global' in prefix_lower:
        score += 20.0
    if preferred_prefix == 'evn' and not _is_global(observation) and '+global' not in prefix_lower:
        score += 5.0

    # Polarization/IF match for LBA-style files.
    if parsed['pols'] is not None and observation.polarizations is not None:
        if parsed['pols'] == observation.polarizations:
            score += 10.0
    if parsed['ifs'] is not None and observation.subbands is not None:
        if parsed['ifs'] == observation.subbands:
            score += 10.0

    return score


def guess_setup_file(observation: Observation) -> Optional[str]:
    """Guess a pySCHED setup file for the observation.

    Looks in ``~/.pysched/setups`` for files whose names match the observing
    band, array and data rate, and returns the best candidate.  Returns None
    if no pySCHED setup directory exists or no suitable file is found.

    Parameters
    ----------
    observation : Observation
        The observation whose parameters drive the search.

    Returns
    -------
    str or None
        Name of the best matching setup file, or None.
    """
    setup_dir = _pysched_setup_dir()
    if not setup_dir.is_dir():
        return None

    candidates = list(setup_dir.glob('*.set'))
    if not candidates:
        return None

    preferred_prefix = _network_prefix(observation)
    scored: list[tuple[float, Path]] = []
    for path in candidates:
        parsed = _parse_setup_filename(path.name)
        score = _setup_file_score(parsed, observation, preferred_prefix)
        if score >= 0.0:
            scored.append((score, path))

    if not scored:
        return None

    scored.sort(key=lambda x: x[0], reverse=True)
    return scored[0][1].name


def _intent_str(stype: SourceType) -> str:
    """Map a SourceType to a SCHED intent string (empty if none applies).

    Parameters
    ----------
    stype : SourceType
        Source type to map.

    Returns
    -------
    str
        SCHED intent string.
    """
    mapping = {
        SourceType.FRINGEFINDER: 'FRINGE_FINDER',
        SourceType.AMPLITUDECAL: 'EMERLIN_AMP',
        SourceType.POLCAL: 'POLCAL',
        SourceType.CHECKSOURCE: 'CHECK',
        SourceType.TARGET: 'TARGET',
        SourceType.PHASECAL: 'PHASE_CAL',
    }
    return mapping.get(stype, '')


# ======================================================================
# Data classes
# ======================================================================

@dataclass
class ScheduledScanBlock:
    """A scan block scheduled at a specific time with timing and quality metrics.

    Attributes
    ----------
    name : str
        Identifier (e.g. 'FF_1', 'R1_D').
    block : ScanBlock
        The original ScanBlock being scheduled.
    start_time : Time
        Start time of this scheduled block.
    end_time : Time
        End time of this scheduled block.
    scans : list[Scan]
        Expanded list of scans filling this time slot.
    n_antennas : int
        Number of antennas that can observe during this block.
    mean_elevation : float
        Mean elevation across observing antennas (degrees).
    """
    name: str
    block: ScanBlock
    start_time: Time
    end_time: Time
    scans: list[Scan] = field(default_factory=list)
    n_antennas: int = 0
    mean_elevation: float = 0.0

    @property
    def duration(self) -> u.Quantity:
        """Duration of this scheduled block."""
        return (self.end_time - self.start_time).to(u.min)


# ======================================================================
# Scheduler
# ======================================================================

class ObservationScheduler:
    """Schedule VLBI scan blocks across the observation time range.

    Parameters
    ----------
    observation : Observation
        The observation containing scan blocks to schedule.
    min_antennas : int
        Minimum antennas required for a valid time slot.
    require_all_antennas : bool
        If True, only schedule when all antennas can observe.
    fringefinder_spec : list[str] or None
        Source names, or ``['N']`` to auto-select *N* FF sources.
    polcal : bool
        Whether polarisation calibration scans are required.
    """

    FF_DUR = 5 * u.min
    FF_INTERVAL = 2 * u.h
    POLCAL_DUR = 5 * u.min
    EMERLIN_3C286_DUR = 5 * u.min
    GAP_DUR = 30 * u.s
    GAP_INTERVAL = 10 * u.min

    # --- CP-SAT objective weights (integer coefficients) ---
    # Spread term rewards a large minimum gap (in grid cells) between
    # consecutive fringe-finder scans.  Separation term penalises the angular
    # distance (deg) between a fringe finder and the science target it
    # interrupts.  Elevation term mildly prefers high-elevation FF placements.
    CP_W_SPREAD = 1000
    CP_W_SEPARATION = 10
    CP_W_ELEVATION = 1
    CP_MAX_SOLVE_SECONDS = 5.0

    def __init__(self, observation: Observation, min_antennas: int = 2,
                 require_all_antennas: bool = False,
                 fringefinder_spec: Optional[list[str]] = None, polcal: bool = False):
        self.obs = observation
        self.min_ant = min_antennas
        self.require_all = require_all_antennas
        self._ff_spec = fringefinder_spec or ['2']
        self._ff_n_scans: Optional[int] = None
        if len(self._ff_spec) == 1 and self._ff_spec[0].isdigit():
            self._ff_n_scans = int(self._ff_spec[0])
        self._polcal = polcal
        self._scheduled: list[ScheduledScanBlock] = []
        # Working copy of the scan blocks: calibrator blocks added while scheduling live only here, so the
        # user's Observation (and its cached source list / visibility / rms) is never mutated.
        self._scans: dict[str, ScanBlock] = dict(observation.scans or {})
        self._rfc_ff_catalog: Optional[calibrators.RFCCatalog] = None
        self._rfc_all_catalog: Optional[calibrators.RFCCatalog] = None
        self._precompute()

    # ------------------------------------------------------------------
    # Pre-computation
    # ------------------------------------------------------------------

    def _precompute(self):
        """Build per-block visibility and mean-elevation arrays over obs.times.

        Precomputes visibility counts and mean elevations for all scan blocks
        at all observation time steps for efficient scheduling.
        """
        self._blocks = list(self._scans.keys())
        self._n_times = len(self.obs.times)
        self._n_ant = len(self.obs.stations)
        self._dt = (self.obs.times[1] - self.obs.times[0]).to(u.min).value
        is_obs = self.obs.is_observable()
        elevs = self.obs.elevations()
        self._vis: dict[str, np.ndarray] = {}
        self._elev: dict[str, np.ndarray] = {}
        self._is_ff: dict[str, bool] = {}
        for name in self._blocks:
            block = self._scans[name]
            self._is_ff[name] = block.has(SourceType.FRINGEFINDER)
            self._vis[name] = (np.sum(np.array(list(is_obs[name].values())), axis=0)
                               if name in is_obs else np.zeros(self._n_times, dtype=int))
            targets = block.sources(SourceType.TARGET) or block.sources(SourceType.FRINGEFINDER) or block.sources()
            if targets and targets[0].name in elevs:
                self._elev[name] = np.nanmean(
                    [elevs[targets[0].name][a].value for a in elevs[targets[0].name]], axis=0)
            else:
                self._elev[name] = np.zeros(self._n_times)

    def _ff_blocks(self) -> list[str]:
        """Get list of fringe-finder block names.

        Returns
        -------
        list[str]
            Names of blocks containing fringe-finder sources.
        """
        return [n for n in self._blocks if self._is_ff.get(n, False)]

    def _sci_blocks(self) -> list[str]:
        """Get list of science block names.

        Returns
        -------
        list[str]
            Names of blocks not containing fringe-finder sources.
        """
        return [n for n in self._blocks if not self._is_ff.get(n, False)]

    def _t2i(self, t: Time) -> int:
        """Convert a Time to the nearest index in obs.times.

        Parameters
        ----------
        t : Time
        Time to convert.

        Returns
        -------
        int
        Index in obs.times array.
        """
        return int(np.argmin(np.abs(self.obs.times.mjd - t.mjd)))

    def _i2t(self, i: int) -> Time:
        """Convert an index to a Time from obs.times.

        Parameters
        ----------
        i : int
        Index in obs.times array.

        Returns
        -------
        Time
        Time at the given index.
        """
        return self.obs.times[min(i, self._n_times - 1)]

    def _vis_at(self, name: str, t: Time) -> int:
        """Get number of visible antennas for a block at a specific time.

        Parameters
        ----------
        name : str
        Block name.
        t : Time
        Observation time.

        Returns
        -------
        int
        Number of visible antennas.
        """
        return int(self._vis[name][self._t2i(t)])

    def _elev_at(self, name: str, t: Time) -> float:
        """Get mean elevation for a block at a specific time.

        Parameters
        ----------
        name : str
        Block name.
        t : Time
        Observation time.

        Returns
        -------
        float
        Mean elevation in degrees.
        """
        return float(self._elev[name][self._t2i(t)])

    def _vis_mean(self, name: str, t0: Time, t1: Time) -> float:
        """Mean number of antennas visible for a block between two times.

        Parameters
        ----------
        name : str
        Block name.
        t0 : Time
        Start time.
        t1 : Time
        End time.

        Returns
        -------
        float
        Mean number of visible antennas.
        """
        i0, i1 = self._t2i(t0), self._t2i(t1)
        arr = self._vis[name][i0:max(i1, i0 + 1)]
        return float(np.mean(arr))

    def _elev_mean(self, name: str, t0: Time, t1: Time) -> float:
        """Mean elevation for a block between two times.

        Parameters
        ----------
        name : str
        Block name.
        t0 : Time
        Start time.
        t1 : Time
        End time.

        Returns
        -------
        float
        Mean elevation in degrees.
        """
        i0, i1 = self._t2i(t0), self._t2i(t1)
        arr = self._elev[name][i0:max(i1, i0 + 1)]
        return float(np.nanmean(arr))

    # ------------------------------------------------------------------
    # Slot helpers
    # ------------------------------------------------------------------

    def _find_best(self, name: str, t0: Time, t1: Time,
                   dur: u.Quantity) -> Optional[tuple[Time, Time, int, float]]:
        """Find optimal placement for a block in a time range with given duration.

        Parameters
        ----------
        name : str
        Block name.
        t0 : Time
        Start of search window.
        t1 : Time
        End of search window.
        dur : Quantity
        Required duration.

        Returns
        -------
        tuple[Time, Time, int, float] or None
        (start, end, min_antennas, mean_elevation) or None if no valid placement.
        """
        i0, i1 = self._t2i(t0), self._t2i(t1)
        steps = max(1, int(np.ceil(dur.to(u.min).value / self._dt)))
        if i1 - i0 < steps:
            return None
        vis, elev = self._vis[name], self._elev[name]
        min_req = self._n_ant if self.require_all else self.min_ant
        # Window k starts at index i0 + k; same min/mean per window as a per-slice loop.
        vis_win = np.lib.stride_tricks.sliding_window_view(vis[i0:i1], steps)
        elev_win = np.lib.stride_tricks.sliding_window_view(elev[i0:i1], steps)
        n_min = vis_win.min(axis=1)
        e_mean = elev_win.mean(axis=1)
        score = n_min * 100.0 + e_mean
        # Invariant: a valid window must also beat the legacy initial best score of -1.0.
        valid = (n_min >= min_req) & ~np.isnan(e_mean) & (score > -1.0)
        if not np.any(valid):
            return None
        best_i = i0 + int(np.argmax(np.where(valid, score, -np.inf)))
        return (self._i2t(best_i), self._i2t(best_i + steps),
                int(np.min(vis[best_i:best_i + steps])),
                float(np.mean(elev[best_i:best_i + steps])))

    def _get_slots(self, sched: list[ScheduledScanBlock], t0: Time, t1: Time) -> list[tuple[Time, Time]]:
        """Return free time slots between already-scheduled blocks.

        Parameters
        ----------
        sched : list[ScheduledScanBlock]
        Already scheduled blocks.
        t0 : Time
        Observation start time.
        t1 : Time
        Observation end time.

        Returns
        -------
        list[tuple[Time, Time]]
        List of (start, end) time tuples for free slots.
        """
        if not sched:
            return [(t0, t1)]
        s = sorted(sched, key=lambda b: b.start_time.mjd)
        slots: list[tuple[Time, Time]] = []
        if s[0].start_time > t0:
            slots.append((t0, s[0].start_time))
        for i in range(len(s) - 1):
            if s[i + 1].start_time > s[i].end_time:
                slots.append((s[i].end_time, s[i + 1].start_time))
        if s[-1].end_time < t1:
            slots.append((s[-1].end_time, t1))
        return slots

    def _register_block(self, name: str):
        """Register a dynamically-added block (FF/polcal/eMERLIN) in the pre-computed arrays.

        Visibility uses the same ERFA batch helper as ``Observation.is_observable`` (a station sees the
        block only when it sees all of its sources); the elevation is the mean over stations of the
        main source elevation, like the targets in ``_precompute``.

        Parameters
        ----------
        name : str
            Block name in ``self._scans`` to register. No-op if already registered.
        """
        if name in self._vis:
            return
        self._blocks.append(name)
        block = self._scans[name]
        self._is_ff[name] = block.has(SourceType.FRINGEFINDER)
        stations = list(self.obs.stations)
        block_sources = block.sources()
        batch_vis = Observation._batch_visibility_erfa(stations, block_sources, self.obs.times)
        n_visible = np.zeros(self._n_times, dtype=int)
        for station in stations:
            station_ok = np.ones(self._n_times, dtype=bool)
            for src in block_sources:
                station_ok &= np.asarray(batch_vis[src.name][station.codename], dtype=bool)
            n_visible += station_ok
        self._vis[name] = n_visible
        main = block.sources(SourceType.TARGET) or block.sources(SourceType.FRINGEFINDER) or block_sources
        if main and stations:
            self._elev[name] = np.nanmean([st.altaz(self.obs.times, main[0]).alt.deg for st in stations], axis=0)
        else:
            self._elev[name] = np.zeros(self._n_times)
        log.debug("Registered block %s: max antennas %d, elevation range %.1f-%.1f deg", name,
                  int(n_visible.max()) if n_visible.size else 0, float(np.nanmin(self._elev[name])),
                  float(np.nanmax(self._elev[name])))

    def _rfc_ff(self) -> calibrators.RFCCatalog:
        """Return the RFC catalog used for automatic fringe-finder selection (>= 0.5 Jy), loaded once.

        Returns
        -------
        calibrators.RFCCatalog
            Catalog with the same filters as the default of ``calibrators.get_fringe_finder_sources``.
        """
        if self._rfc_ff_catalog is None:
            self._rfc_ff_catalog = calibrators.RFCCatalog(min_flux=0.5 * u.Jy, band='c')
        return self._rfc_ff_catalog

    def _rfc_all(self) -> calibrators.RFCCatalog:
        """Return the unfiltered RFC catalog (min flux 0 Jy) used for name lookups, loaded once.

        Returns
        -------
        calibrators.RFCCatalog
            Catalog used to resolve named fringe finders and polcal sources.
        """
        if self._rfc_all_catalog is None:
            self._rfc_all_catalog = calibrators.RFCCatalog(min_flux=0.0 * u.Jy, band='c')
        return self._rfc_all_catalog

    # ------------------------------------------------------------------
    # Fringe-finder selection  (NEW)
    # ------------------------------------------------------------------

    def _select_fringefinders(self) -> list[str]:
        """Select FF sources and register them as scan blocks.

        Strategy
        --------
        1. If the user gave explicit source names, use those.
        2. Otherwise find FFs visible by ALL antennas for the ENTIRE observation;
           pick the one closest to the mean science-target position.
        3. Fallback: find FFs visible by at least *min_ant* antennas; for each
           FF time slot pick the best-visible source closest to the targets.

        Returns
        -------
        list[str]
            Block names registered in ``self._scans``.
        """
        existing = self._ff_blocks()
        if existing:
            return existing

        spec = self._ff_spec
        if spec and not (len(spec) == 1 and spec[0].isdigit()):
            return self._register_ff_sources(self._lookup_ff_sources(spec))

        # Try all-antenna, all-time visibility first
        cands, _, _ = calibrators.get_fringe_finder_sources(
            self.obs.stations, self.obs.times, min_elevation=20 * u.deg, min_flux=0.5 * u.Jy,
            catalog=self._rfc_ff(), require_all_stations=True)

        if cands:
            best = self._pick_closest_ff(cands)
            return self._register_ff_sources([best])

        # Relax to partial visibility
        cands, _, _ = calibrators.get_fringe_finder_sources(
            self.obs.stations, self.obs.times, min_elevation=20 * u.deg, min_flux=0.5 * u.Jy,
            catalog=self._rfc_ff(), require_all_stations=False)

        if cands:
            best = self._pick_closest_ff(cands)
            return self._register_ff_sources([best])

        log.warning("No fringe-finder sources found for this observation.")
        return []

    def _pick_closest_ff(self, candidates: list) -> Source:
        """From a list of CalibratorSource candidates, pick the one closest to science targets.

        Parameters
        ----------
        candidates : list
        List of CalibratorSource candidates.

        Returns
        -------
        Source
        The closest source converted to a Source object.
        """
        mean_coord = self._mean_target_coord()
        if mean_coord is None:
            src = candidates[0]
            return Source(name=src.name, coordinates=src.coord,
                          source_type=SourceType.FRINGEFINDER, other_names=[src.ivsname])
        best_src, best_sep = candidates[0], float('inf')
        for c in candidates:
            sep = c.coord.separation(mean_coord).deg
            if sep < best_sep:
                best_sep, best_src = sep, c
        return Source(name=best_src.name, coordinates=best_src.coord,
                      source_type=SourceType.FRINGEFINDER, other_names=[best_src.ivsname])

    def _mean_target_coord(self) -> Optional[SkyCoord]:
        """Return the mean sky position of all science targets, or None.

        Returns
        -------
        SkyCoord or None
        Mean coordinate of all science targets, or None if no targets.
        """
        ras, decs = [], []
        for name in self._sci_blocks():
            for src in self._scans[name].sources(SourceType.TARGET):
                ras.append(src.coord.ra.deg)
                decs.append(src.coord.dec.deg)
        if not ras:
            return None
        return SkyCoord(float(np.mean(ras)), float(np.mean(decs)), unit=u.deg)

    def _register_ff_sources(self, sources: list[Source]) -> list[str]:
        """Register unique FF Source objects as scan blocks.

        Parameters
        ----------
        sources : list[Source]
        List of Source objects to register.

        Returns
        -------
        list[str]
        Block names registered in the scheduler working blocks (``self._scans``).
        """
        seen: set[str] = set()
        names: list[str] = []
        for src in sources:
            if src.name in seen:
                continue
            seen.add(src.name)
            block_name = f"FF_{src.name}"
            if block_name not in self._scans:
                self._scans[block_name] = ScanBlock([Scan(src, duration=self.FF_DUR)])
                self._register_block(block_name)
            names.append(block_name)
        return names

    def _lookup_ff_sources(self, names: list[str]) -> list[Source]:
        """Look up named fringe-finder sources from the RFC catalog.

        Each entry in *names* can be 'name/coordinates' to supply explicit
        coordinates, plain coordinates, or a catalog name.

        Parameters
        ----------
        names : list[str]
        List of source specifications.

        Returns
        -------
        list[Source]
        List of Source objects with FRINGEFINDER type.
        """
        cat = self._rfc_all()
        result: list[Source] = []
        for spec in names:
            parsed_name, parsed_coord = Source.parse_source_spec(spec)
            if parsed_name is not None and parsed_coord is not None:
                result.append(Source(
                    parsed_name, coordinates=Source._parse_coord_str(parsed_coord),
                    source_type=SourceType.FRINGEFINDER))
            else:
                lookup = parsed_name or spec
                rfc_src = cat.get_source(lookup)
                if rfc_src is not None:
                    result.append(Source(name=rfc_src.name, coordinates=rfc_src.coord,
                                        source_type=SourceType.FRINGEFINDER,
                                        other_names=[rfc_src.ivsname]))
                else:
                    try:
                        result.append(Source.source_from_str(spec, source_type=SourceType.FRINGEFINDER))
                    except ValueError:
                        log.warning("Fringe-finder source '%s' not found — skipping.", spec)
        return result

    # ------------------------------------------------------------------
    # Polcal / eMERLIN helpers
    # ------------------------------------------------------------------

    def _create_polcal_blocks(self) -> list[str]:
        """Create 3C84 / OQ208 / DA193 polcal blocks.

        Returns
        -------
        list[str]
            Block names registered in the scheduler working blocks (``self._scans``).
        """
        cat = self._rfc_all()
        added: list[str] = []
        for polcal_name in _POLCAL_NAMES:
            rfc_src = cat.get_source(polcal_name)
            if rfc_src is None:
                continue
            src = Source(name=rfc_src.name, coordinates=rfc_src.coord,
                         source_type=SourceType.POLCAL, other_names=[rfc_src.ivsname])
            block_name = f"POLCAL_{rfc_src.ivsname}"
            self._scans[block_name] = ScanBlock([Scan(src, duration=self.POLCAL_DUR)])
            self._register_block(block_name)
            added.append(block_name)
        return added

    def _has_emerlin(self) -> bool:
        """Check if eMERLIN stations (other than Jb2, which is a regular EVN station) are in the array.

        Returns
        -------
        bool
        True if any station in `_EMERLIN_CODES` is present.
        """
        return any(s.codename.upper() in _EMERLIN_CODES for s in self.obs.stations)

    def _has_jb1(self) -> bool:
        """Return True if Jodrell Bank Mk2 (Jb1 / JB1) is in the array.

        Returns
        -------
        bool
        True if JB1 is present.
        """
        return any(s.codename.upper() == 'JB1' for s in self.obs.stations)

    def _create_emerlin_3c286_block(self) -> Optional[str]:
        """Create a 3C286 scan block for eMERLIN flux-scale calibration.

        Returns
        -------
        str or None
        Block name if created, None otherwise.
        """
        src = Source(name='3C286', coordinates=_3C286_COORD, source_type=SourceType.AMPLITUDECAL)
        block_name = 'eMERLIN_3C286'
        self._scans[block_name] = ScanBlock([Scan(src, duration=self.EMERLIN_3C286_DUR)])
        self._register_block(block_name)
        return block_name

    # ------------------------------------------------------------------
    # FF scheduling  (REWRITTEN — evenly spread, never consecutive)
    # ------------------------------------------------------------------

    def _effective_obs_window(self, t0: Time, t1: Time) -> tuple[Time, Time]:
        """Return the effective observation window where science is possible.

        Shrinks [t0, t1] to the range where at least one science block has
        >= min_ant antennas visible.  FF scans should be placed only within
        this window so they bracket actual science time.

        Parameters
        ----------
        t0 : Time
        Observation start time.
        t1 : Time
        Observation end time.

        Returns
        -------
        tuple[Time, Time]
        (effective_start, effective_end).
        """
        sci_names = self._sci_blocks()
        if not sci_names:
            return t0, t1
        # Combine visibility across all science blocks (take max per time step)
        combined = np.zeros(self._n_times, dtype=int)
        for n in sci_names:
            combined = np.maximum(combined, self._vis[n])
        ok = np.where(combined >= self.min_ant)[0]
        if len(ok) == 0:
            return t0, t1
        eff_t0 = self._i2t(int(ok[0]))
        eff_t1 = self._i2t(int(ok[-1]))
        # Pad the end by FF_DUR so we can fit a closing FF scan
        eff_t1 = min(t1, eff_t1 + self.FF_DUR)
        return max(t0, eff_t0), eff_t1

    def _schedule_ff(self, ff_names: list[str], t0: Time, t1: Time) -> list[ScheduledScanBlock]:
        """Place FF scans evenly across the observation, never adjacent.

        One scan at start, one at end, additional every ~FF_INTERVAL in between.
        FF scans are only placed within the effective observation window where
        science targets are observable.

        Parameters
        ----------
        ff_names : list[str]
            Block names for FF sources (typically one element).
        t0 : Time
            Observation start time.
        t1 : Time
            Observation end time.

        Returns
        -------
        list[ScheduledScanBlock]
            Scheduled FF blocks in chronological order.
        """
        if not ff_names:
            return []

        # Restrict to where science is actually observable
        eff_t0, eff_t1 = self._effective_obs_window(t0, t1)
        dur = self.FF_DUR
        obs_dur_h = (eff_t1 - eff_t0).to(u.h).value

        if self._ff_n_scans is not None:
            n_ff = max(1, self._ff_n_scans)
            log.info("Using user-requested FF scan count: %d", n_ff)
        elif obs_dur_h <= 1.5:
            n_ff = 1
        elif obs_dur_h <= 3.0:
            n_ff = 2
        else:
            n_mid = max(0, int((obs_dur_h - dur.to(u.h).value) / self.FF_INTERVAL.to(u.h).value))
            n_ff = n_mid + 1

        total_ff_time = n_ff * dur
        total_sci_time = (eff_t1 - eff_t0) - total_ff_time
        if total_sci_time.to(u.min).value < 0:
            total_sci_time = 0 * u.min

        gap = total_sci_time / max(n_ff - 1, 1) if n_ff > 1 else 0 * u.min
        starts: list[Time] = []
        t = eff_t0
        for _ in range(n_ff):
            starts.append(t)
            t = t + dur + gap

        result: list[ScheduledScanBlock] = []
        for idx, ts in enumerate(starts):
            te = ts + dur
            block_name = ff_names[idx % len(ff_names)]
            n_vis = self._vis_at(block_name, ts)
            if n_vis < self.min_ant:
                log.info("Skipping FF slot at %s — only %d antennas can observe.", ts.iso, n_vis)
                continue
            block = self._scans[block_name]
            result.append(ScheduledScanBlock(
                name=f"FF_{len(result) + 1}", block=block, start_time=ts, end_time=te,
                scans=list(block.scans), n_antennas=n_vis,
                mean_elevation=self._elev_at(block_name, ts)))
        return result

    # ------------------------------------------------------------------
    # Polcal / single-block scheduling  (unchanged logic)
    # ------------------------------------------------------------------

    def _schedule_polcal(self, polcal_names: list[str],
                         slots: list[tuple[Time, Time]]) -> list[ScheduledScanBlock]:
        """Spread polcal scans at 10 %, 50 %, 90 % of available time.

        Parameters
        ----------
        polcal_names : list[str]
            Block names for polcal sources.
        slots : list[tuple[Time, Time]]
            Available time slots.

        Returns
        -------
        list[ScheduledScanBlock]
            Scheduled polcal blocks.
        """
        if not polcal_names or not slots:
            return []
        n_polcal = min(3, len(polcal_names))
        dur = self.POLCAL_DUR
        total = sum(((s1 - s0).to(u.min).value for s0, s1 in slots)) * u.min
        if total < dur:
            return []
        placed: list[ScheduledScanBlock] = []
        for frac, pc_name in zip([0.1, 0.5, 0.9][:n_polcal], polcal_names):
            block = self._scans[pc_name]
            elapsed = 0.0 * u.min
            target_time = frac * total
            for s0, s1 in slots:
                slot_dur = (s1 - s0).to(u.min)
                if elapsed + slot_dur >= target_time and slot_dur >= dur:
                    offset = max(target_time - elapsed, 0.0 * u.min)
                    ts = s0 + offset
                    if ts + dur <= s1:
                        placed.append(ScheduledScanBlock(
                            name=f"POLCAL_{pc_name.split('_')[-1]}", block=block,
                            start_time=ts, end_time=ts + dur, scans=[block.scans[0]],
                            n_antennas=self._vis_at(pc_name, ts),
                            mean_elevation=self._elev_at(pc_name, ts)))
                        break
                elapsed += slot_dur
        return placed

    def _schedule_single_block(self, block_name: str, slots: list[tuple[Time, Time]],
                               label: str) -> Optional[ScheduledScanBlock]:
        """Place a single auxiliary block in the best available slot.

        Parameters
        ----------
        block_name : str
        Block name to schedule.
        slots : list[tuple[Time, Time]]
        Available time slots.
        label : str
        Label for the scheduled block.

        Returns
        -------
        ScheduledScanBlock or None
        Scheduled block or None if no valid slot found.
        """
        block = self._scans[block_name]
        dur = sum(s.duration.to(u.min).value for s in block.scans) * u.min
        best = None
        for s0, s1 in slots:
            if (s1 - s0).to(u.min) >= dur:
                r = self._find_best(block_name, s0, s1, dur)
                if r and (not best or r[2] * 100 + r[3] > best[2] * 100 + best[3]):
                    best = r
        if best:
            ts, te, n, e = best
            return ScheduledScanBlock(name=label, block=block, start_time=ts, end_time=te,
                                      scans=[block.scans[0]], n_antennas=n, mean_elevation=e)
        return None

    # ------------------------------------------------------------------
    # Science scheduling  (REWRITTEN — optimise per slot)
    # ------------------------------------------------------------------

    def _schedule_sci(self, names: list[str], slots: list[tuple[Time, Time]]) -> list[ScheduledScanBlock]:
        """Fill each free slot with the science block that has the best observing conditions there.

        For a single science block all slots are filled with it.  For multiple
        blocks the one with the highest (n_antennas * 100 + mean_elevation)
        score wins each slot.

        Parameters
        ----------
        names : list[str]
            Science block names.
        slots : list[tuple[Time, Time]]
            Available time windows between FF / calibrator scans.

        Returns
        -------
        list[ScheduledScanBlock]
            Scheduled science blocks.
        """
        if not names or not slots:
            return []
        scheduled: list[ScheduledScanBlock] = []
        for s0, s1 in slots:
            slot_dur = (s1 - s0).to(u.min)
            best_name: Optional[str] = None
            best_score = -1.0
            for name in names:
                block = self._scans[name]
                min_dur = sum(s.duration.to(u.min).value for s in block.scans) * u.min
                if slot_dur < min_dur:
                    continue
                score = self._vis_mean(name, s0, s1) * 100.0 + self._elev_mean(name, s0, s1)
                if score > best_score:
                    best_score, best_name = score, name
            if best_name is None:
                continue
            block = self._scans[best_name]
            try:
                scans = block.fill(slot_dur)
            except ValueError:
                scans = list(block.scans)
            mid = s0 + (s1 - s0) * 0.5
            scheduled.append(ScheduledScanBlock(
                name=best_name, block=block, start_time=s0, end_time=s1, scans=scans,
                n_antennas=self._vis_at(best_name, mid),
                mean_elevation=self._elev_at(best_name, mid)))
        return scheduled

    # ------------------------------------------------------------------
    # Multi-source scheduling helpers
    # ------------------------------------------------------------------

    def _source_main_coord(self, name: str) -> Optional[SkyCoord]:
        """Return the sky coordinate of the primary science target in a scan block.

        Parameters
        ----------
        name : str
            Block name in ``self._scans``.

        Returns
        -------
        SkyCoord or None
            Coordinate of the primary target, or None if no target found.
        """
        block = self._scans.get(name)
        if block is None:
            return None
        targets = block.sources(SourceType.TARGET) or block.sources(SourceType.PULSAR) or block.sources()
        return targets[0].coord if targets else None

    def _find_optimal_window(self, name: str, t0: Time, t1: Time, dur: u.Quantity) -> tuple[Time, float]:
        """Find the start time of the best observing window for a block.

        Searches the interval [t0, t1] for the placement of *dur* that maximises
        n_antennas * 100 + mean_elevation.

        Parameters
        ----------
        name : str
        Block name.
        t0 : Time
        Start of search window.
        t1 : Time
        End of search window.
        dur : Quantity
        Required duration.

        Returns
        -------
        tuple[Time, float]
        (best_start, score) — score is -1 if no valid window found.
        """
        result = self._find_best(name, t0, t1, dur)
        if result:
            ts, _te, n, e = result
            return ts, float(n) * 100.0 + e
        return t0, -1.0

    def _nearest_neighbor_order(self, names: list[str]) -> list[str]:
        """Reorder names using nearest-neighbour heuristic to minimise total angular slewing.

        Starts from the first element of *names* (which should already be sorted
        by optimal observation window) and greedily picks the geographically
        closest remaining source at each step.

        Parameters
        ----------
        names : list[str]
            Source block names, pre-sorted by optimal start time.

        Returns
        -------
        list[str]
            Reordered names.
        """
        if len(names) <= 2:
            return list(names)
        coords: dict[str, Optional[SkyCoord]] = {n: self._source_main_coord(n) for n in names}
        remaining = list(names)
        ordered = [remaining.pop(0)]
        while remaining:
            last_coord = coords.get(ordered[-1])
            if last_coord is None:
                ordered.append(remaining.pop(0))
                continue
            best_next = min(remaining, key=lambda n: (
                last_coord.separation(coords[n]).deg if coords.get(n) is not None else 360.0))
            ordered.append(best_next)
            remaining.remove(best_next)
        return ordered

    def _choose_ff_positions(self, n_sources: int, n_ff: int) -> set[int]:
        """Choose at which source boundaries to insert FF scans.

        Positions range from 0 (before source 0) to n_sources (after last source).
        FFs are spread as evenly as possible across the observation.

        Parameters
        ----------
        n_sources : int
        Number of science sources.
        n_ff : int
        Number of FF scans to place.

        Returns
        -------
        set[int]
        Set of boundary positions where FF scans should be inserted.
        """
        if n_ff <= 0:
            return set()
        if n_ff >= n_sources + 1:
            return set(range(n_sources + 1))
        positions: set[int] = set()
        for k in range(n_ff):
            idx = round((k + 0.5) * (n_sources + 1) / n_ff)
            positions.add(min(max(idx, 0), n_sources))
        return positions

    def _determine_ff_count(self, obs_dur_h: float) -> int:
        """Determine the number of FF scans for a given observation duration.

        Respects user-specified count.  Otherwise:
        <= 1.5 h → 1, <= 3 h → 2, > 3 h → 2 + one per FF_INTERVAL.

        Parameters
        ----------
        obs_dur_h : float
            Observation duration in hours.

        Returns
        -------
        int
            Number of FF scans to schedule.
        """
        if self._ff_n_scans is not None:
            return max(1, self._ff_n_scans)
        if obs_dur_h <= 1.5:
            return 1
        if obs_dur_h <= 3.0:
            return 2
        extra = max(0, int(obs_dur_h / self.FF_INTERVAL.to(u.h).value) - 1)
        return 2 + extra

    def _schedule_multi_source(self, sci_names: list[str], ff_names: list[str],
                                t0: Time, t1: Time) -> list[ScheduledScanBlock]:
        """Schedule multiple science sources with equalised time and interleaved FF scans.

        Algorithm
        ---------
        1. Calculate equal science time per source (total_time - ff_overhead).
        2. Find the optimal observing window for each source.
        3. Sort sources chronologically by their optimal window; break ties with
           nearest-neighbour sky distance to minimise slewing.
        4. Place FF scans at evenly-distributed boundaries between source blocks.
        5. Place each source block sequentially; the duration is adjusted so that
           the last block ends exactly at t1.

        Parameters
        ----------
        sci_names : list[str]
            Science block names (already filtered, no FF/polcal/eMERLIN blocks).
        ff_names : list[str]
            FF block names registered in ``self._scans``.
        t0 : Time
            Observation start time.
        t1 : Time
            Observation end time.

        Returns
        -------
        list[ScheduledScanBlock]
            All scheduled blocks (science + FF) in chronological order.
        """
        n_sci = len(sci_names)
        obs_dur_h = (t1 - t0).to(u.h).value
        n_ff = self._determine_ff_count(obs_dur_h)

        # --- time budget ---
        total_ff_time = n_ff * self.FF_DUR
        sci_time_total = max((t1 - t0) - total_ff_time, (t1 - t0) * 0.1)
        time_per_source = sci_time_total / n_sci

        # --- find optimal windows ---
        source_peak: dict[str, Time] = {}
        for name in sci_names:
            peak_start, _ = self._find_optimal_window(name, t0, t1, time_per_source)
            source_peak[name] = peak_start

        # --- sort by optimal window start, then minimise slewing ---
        ordered = sorted(sci_names, key=lambda n: source_peak[n].mjd)
        ordered = self._nearest_neighbor_order(ordered)

        # --- choose where to insert FF scans ---
        ff_positions = self._choose_ff_positions(n_sci, n_ff)
        log.info("Multi-source schedule: %d sources, %d FF scans at positions %s",
                 n_sci, n_ff, sorted(ff_positions))

        # --- build timeline ---
        all_blocks: list[ScheduledScanBlock] = []
        ff_count = 0
        current_t = t0

        def _place_ff(t: Time, count: int) -> Time:
            """Place one FF block starting at t; returns new current time."""
            if not ff_names:
                return t
            fn = ff_names[count % len(ff_names)]
            ff_block = self._scans[fn]
            n_vis = self._vis_at(fn, t)
            all_blocks.append(ScheduledScanBlock(
                name=f"FF_{count + 1}", block=ff_block,
                start_time=t, end_time=t + self.FF_DUR,
                scans=list(ff_block.scans), n_antennas=n_vis,
                mean_elevation=self._elev_at(fn, t)))
            return t + self.FF_DUR

        for i, name in enumerate(ordered):
            # Possibly insert FF before this source
            if i in ff_positions and ff_count < n_ff:
                current_t = _place_ff(current_t, ff_count)
                ff_count += 1

            # Allocate time: divide remaining time equally among remaining sources,
            # reserving time for any FFs still to be placed after this source.
            remaining_sci_blocks = n_sci - i
            remaining_ff_time = (n_ff - ff_count) * self.FF_DUR
            available_to_end = (t1 - current_t) - remaining_ff_time
            this_dur = max(available_to_end / remaining_sci_blocks, 1.0 * u.min)

            block = self._scans[name]
            sci_dur = min(this_dur, (t1 - current_t)).to(u.min)
            try:
                scans = block.fill(sci_dur)
            except (ValueError, Exception):
                scans = list(block.scans)

            # Use actual scan durations for precise time tracking (fill() may not use the full sci_dur)
            actual_dur = sum(s.duration.to(u.min).value for s in scans) * u.min
            sci_end = current_t + actual_dur
            mid = current_t + actual_dur * 0.5
            all_blocks.append(ScheduledScanBlock(
                name=name, block=block, start_time=current_t, end_time=sci_end,
                scans=scans, n_antennas=self._vis_at(name, mid),
                mean_elevation=self._elev_at(name, mid)))
            current_t = sci_end

        # Place any remaining FFs at the end
        while ff_count < n_ff and ff_names and current_t + self.FF_DUR <= t1:
            current_t = _place_ff(current_t, ff_count)
            ff_count += 1

        return sorted(all_blocks, key=lambda b: b.start_time.mjd)

    # ------------------------------------------------------------------
    # CP-SAT scheduling  (constraint-programming layout)
    # ------------------------------------------------------------------

    def _block_source_coord(self, name: str) -> Optional[SkyCoord]:
        """Return the sky coordinate of the first source in a scan block.

        Parameters
        ----------
        name : str
        Block name.

        Returns
        -------
        SkyCoord or None
        Coordinate of the first source, or None if block not found.
        """
        block = self._scans.get(name)
        if block is None:
            return None
        srcs = block.sources()
        return srcs[0].coord if srcs else None

    def _min_block_duration(self, name: str) -> u.Quantity:
        """Return the minimum duration to fit one pass of a scan block.

        Parameters
        ----------
        name : str
        Block name.

        Returns
        -------
        Quantity
        Minimum duration in minutes.
        """
        block = self._scans[name]
        return sum(s.duration.to(u.min).value for s in block.scans) * u.min

    def _layout_science(self, sci_names: list[str], eff_t0: Time, eff_t1: Time) -> list[ScheduledScanBlock]:
        """Lay out science blocks contiguously across [eff_t0, eff_t1] (no calibrators).

        For a single science block the whole window is filled with it.  For
        multiple blocks the time is equalised per source, ordered by optimal
        observing window, and slew-minimised via nearest-neighbour.  Fringe
        finders and other calibrators are inserted afterwards by splitting these
        science blocks (see ``_insert_aux``).

        Parameters
        ----------
        sci_names : list[str]
            Science block names.
        eff_t0 : Time
            Effective observation window start where science is observable.
        eff_t1 : Time
            Effective observation window end where science is observable.

        Returns
        -------
        list[ScheduledScanBlock]
            Contiguous science blocks covering the window.
        """
        if not sci_names:
            return []
        if len(sci_names) == 1:
            name = sci_names[0]
            block = self._scans[name]
            span = (eff_t1 - eff_t0)
            try:
                scans = block.fill(span.to(u.min))
            except Exception:
                scans = list(block.scans)
            mid = eff_t0 + span * 0.5
            return [ScheduledScanBlock(
                name=name, block=block, start_time=eff_t0, end_time=eff_t1, scans=scans,
                n_antennas=self._vis_at(name, mid), mean_elevation=self._elev_at(name, mid))]

        n = len(sci_names)
        tps = (eff_t1 - eff_t0) / n
        peak = {nm: self._find_optimal_window(nm, eff_t0, eff_t1, tps)[0] for nm in sci_names}
        ordered = sorted(sci_names, key=lambda nm: peak[nm].mjd)
        ordered = self._nearest_neighbor_order(ordered)

        blocks: list[ScheduledScanBlock] = []
        cur = eff_t0
        for i, name in enumerate(ordered):
            remaining = n - i
            this_dur = (eff_t1 - cur) / remaining
            block = self._scans[name]
            sci_dur = min(this_dur, (eff_t1 - cur)).to(u.min)
            try:
                scans = block.fill(sci_dur)
            except Exception:
                scans = list(block.scans)
            actual = sum(s.duration.to(u.min).value for s in scans) * u.min
            if actual.to(u.min).value <= 0:
                actual = sci_dur
            end = cur + actual
            mid = cur + actual * 0.5
            blocks.append(ScheduledScanBlock(
                name=name, block=block, start_time=cur, end_time=end, scans=scans,
                n_antennas=self._vis_at(name, mid), mean_elevation=self._elev_at(name, mid)))
            cur = end
        return blocks

    def _occupant_at(self, layout: list[ScheduledScanBlock], t: Time) -> Optional[tuple[str, Optional[SkyCoord]]]:
        """Return (block_name, target_coord) of the science block covering time t, or None.

        Parameters
        ----------
        layout : list[ScheduledScanBlock]
        Scheduled science blocks.
        t : Time
        Time to check.

        Returns
        -------
        tuple[str, SkyCoord] or None
        (block_name, target_coord) if t is within a block, else None.
        """
        for sb in layout:
            if sb.start_time <= t < sb.end_time:
                return sb.name, self._source_main_coord(sb.name)
        return None

    def _cpsat_ff_placement(self, ff_names: list[str], layout: list[ScheduledScanBlock],
                            eff_t0: Time, eff_t1: Time) -> list[tuple[Time, str]]:
        """Decide fringe-finder times and sources via CP-SAT optimisation.

        The observation window is discretised into ``FF_DUR``-sized cells.  Each
        FF scan is assigned to a cell (and a source) such that:

        - exactly ``n_ff`` scans are placed (from ``_determine_ff_count``);
        - each FF lies fully inside one science block (so it cleanly splits it);
        - the FF source is visible by the required number of antennas there;
        - the **minimum gap between consecutive FF scans is maximised** (even
          spread along the observation, ``CP_W_SPREAD``);
        - the **angular separation between each FF and the science target it
          interrupts is minimised** (``CP_W_SEPARATION``);
        - high FF elevation is mildly preferred (``CP_W_ELEVATION``).

        Parameters
        ----------
        ff_names : list[str]
            Candidate FF block names (registered in ``self._scans``).
        layout : list[ScheduledScanBlock]
            The science-only layout (provides the interrupted-target coords).
        eff_t0 : Time
            Effective observation window start.
        eff_t1 : Time
            Effective observation window end.

        Returns
        -------
        list[tuple[Time, str]]
            (start_time, ff_block_name) pairs in chronological order.  Empty if
            ortools is unavailable or no feasible placement exists.
        """
        if not _HAS_ORTOOLS or not ff_names or not layout:
            return []
        dur = self.FF_DUR
        cell_min = dur.to(u.min).value
        win_min = (eff_t1 - eff_t0).to(u.min).value
        K = int(win_min // cell_min)
        if K < 1:
            return []

        cell_t = [eff_t0 + (c * cell_min) * u.min for c in range(K)]
        min_req = self._n_ant if self.require_all else self.min_ant
        ff_coord = {fn: self._block_source_coord(fn) for fn in ff_names}

        # Candidate (cell, ff-source) placements with separation + elevation cost.
        cand: dict[tuple[int, int], tuple[float, float]] = {}
        cand_by_cell: dict[int, list[int]] = {}
        sep_cache: dict[tuple[str, str], float] = {}
        for c in range(K):
            occ_s = self._occupant_at(layout, cell_t[c])
            occ_e = self._occupant_at(layout, cell_t[c] + dur - 1 * u.s)
            if occ_s is None or occ_e is None or occ_s[0] != occ_e[0]:
                continue  # FF would cross a science-block boundary; skip cell
            occ_coord = occ_s[1]
            for fi, fn in enumerate(ff_names):
                if self._vis_at(fn, cell_t[c]) < min_req:
                    continue
                sep = 0.0
                if occ_coord is not None and ff_coord[fn] is not None:
                    # The occupant coordinate depends only on its block name, so cache per (FF, occupant).
                    sep_key = (fn, occ_s[0])
                    if sep_key not in sep_cache:
                        sep_cache[sep_key] = float(ff_coord[fn].separation(occ_coord).deg)
                    sep = sep_cache[sep_key]
                cand[(c, fi)] = (sep, self._elev_at(fn, cell_t[c]))
                cand_by_cell.setdefault(c, []).append(fi)
        if not cand:
            return []

        usable_cells = sorted(cand_by_cell.keys())
        n_ff = min(self._determine_ff_count((eff_t1 - eff_t0).to(u.h).value), len(usable_cells))
        if n_ff < 1:
            return []

        model = cp_model.CpModel()
        x: dict[tuple[int, int, int], object] = {}
        for i in range(n_ff):
            for (c, fi) in cand:
                x[i, c, fi] = model.NewBoolVar(f"x_{i}_{c}_{fi}")
            model.Add(sum(x[i, c, fi] for (c, fi) in cand) == 1)
        for c in usable_cells:
            model.Add(sum(x[i, c, fi] for i in range(n_ff) for fi in cand_by_cell[c]) <= 1)

        pos = [model.NewIntVar(0, K - 1, f"pos_{i}") for i in range(n_ff)]
        for i in range(n_ff):
            model.Add(pos[i] == sum(c * x[i, c, fi] for (c, fi) in cand))
        for i in range(n_ff - 1):
            model.Add(pos[i + 1] >= pos[i] + 1)

        obj = []
        if n_ff > 1:
            gmin = model.NewIntVar(0, K, "gmin")
            for i in range(n_ff - 1):
                model.Add(gmin <= pos[i + 1] - pos[i])
            obj.append(self.CP_W_SPREAD * gmin)
        else:
            dev = model.NewIntVar(0, K, "dev")
            model.AddAbsEquality(dev, pos[0] - K // 2)
            obj.append(-self.CP_W_SPREAD * dev)

        sep_term = sum(int(round(cand[(c, fi)][0])) * x[i, c, fi]
                       for i in range(n_ff) for (c, fi) in cand)
        elev_term = sum(int(round(cand[(c, fi)][1])) * x[i, c, fi]
                        for i in range(n_ff) for (c, fi) in cand)
        model.Maximize(sum(obj) - self.CP_W_SEPARATION * sep_term + self.CP_W_ELEVATION * elev_term)

        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = self.CP_MAX_SOLVE_SECONDS
        solver.parameters.random_seed = 0
        solver.parameters.num_search_workers = 1
        status = solver.Solve(model)
        if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            log.warning("CP-SAT FF placement found no solution (status %s).", solver.StatusName(status))
            return []

        placements: list[tuple[Time, str]] = []
        for i in range(n_ff):
            for (c, fi) in cand:
                if solver.Value(x[i, c, fi]) == 1:
                    placements.append((cell_t[c], ff_names[fi]))
                    break
        placements.sort(key=lambda p: p[0].mjd)
        log.info("CP-SAT placed %d FF scans (spread=%s, status=%s).",
                 len(placements), 'maximin' if n_ff > 1 else 'centred', solver.StatusName(status))
        return placements

    def _polcal_target_times(self, eff_t0: Time, eff_t1: Time) -> list[tuple[Time, str]]:
        """Create polcal blocks and return their (time, name) placements near 10/50/90% of the window.

        Each polcal is placed at the grid time closest to its target fraction where the polcal source
        is visible by the required number of antennas; polcals never visible in the window are skipped.

        Parameters
        ----------
        eff_t0 : Time
        Effective observation window start.
        eff_t1 : Time
        Effective observation window end.

        Returns
        -------
        list[tuple[Time, str]]
        List of (start_time, block_name) tuples.
        """
        names = self._create_polcal_blocks()
        total = (eff_t1 - eff_t0)
        min_req = self._n_ant if self.require_all else self.min_ant
        i0, i_last = self._t2i(eff_t0), self._t2i(eff_t1 - self.POLCAL_DUR)
        result: list[tuple[Time, str]] = []
        for frac, pc in zip([0.1, 0.5, 0.9], names):
            target = eff_t0 + total * frac
            if self._vis_at(pc, target) >= min_req:
                result.append((target, pc))
                continue
            ok_idx = np.nonzero(self._vis[pc][i0:i_last + 1] >= min_req)[0] + i0
            if ok_idx.size == 0:
                log.warning("Polcal block %s is not visible by %d antennas in the window; skipped.", pc, min_req)
                continue
            i_target = self._t2i(target)
            best = int(ok_idx[np.argmin(np.abs(ok_idx - i_target))])
            log.info("Polcal block %s moved from %s to %s (source not visible at target time).",
                     pc, target.iso, self._i2t(best).iso)
            result.append((self._i2t(best), pc))
        return result

    def _make_aux_block(self, name: str, a: Time, b: Time, label: str) -> ScheduledScanBlock:
        """Build a single-scan calibrator ScheduledScanBlock (FF/polcal/eMERLIN).

        Parameters
        ----------
        name : str
        Block name in ``self._scans``.
        a : Time
        Start time.
        b : Time
        End time.
        label : str
        Label for the scheduled block.

        Returns
        -------
        ScheduledScanBlock
        The scheduled calibrator block.
        """
        block = self._scans[name]
        return ScheduledScanBlock(
            name=label, block=block, start_time=a, end_time=b, scans=list(block.scans),
            n_antennas=self._vis_at(name, a), mean_elevation=self._elev_at(name, a))

    def _make_sci_chunk(self, name: str, a: Time, b: Time) -> Optional[ScheduledScanBlock]:
        """Build a science ScheduledScanBlock for [a, b], or None if too short to fit one pass.

        Parameters
        ----------
        name : str
        Block name.
        a : Time
        Start time.
        b : Time
        End time.

        Returns
        -------
        ScheduledScanBlock or None
        Scheduled science block, or None if duration insufficient.
        """
        span = (b - a)
        if span.to(u.min).value <= 0 or span < self._min_block_duration(name):
            return None
        block = self._scans[name]
        try:
            scans = block.fill(span.to(u.min))
        except Exception:
            scans = list(block.scans)
        mid = a + span * 0.5
        return ScheduledScanBlock(
            name=name, block=block, start_time=a, end_time=b, scans=scans,
            n_antennas=self._vis_at(name, mid), mean_elevation=self._elev_at(name, mid))

    def _insert_aux(self, layout: list[ScheduledScanBlock],
                    placements: list[tuple[Time, str, u.Quantity, str]]) -> list[ScheduledScanBlock]:
        """Split science blocks to insert calibrator scans at the requested times.

        Each placement ``(start, block_name, duration, label)`` is assigned to the
        science block that contains its start time and inserted there, splitting
        the surrounding science into the chunks before and after it.  Placements
        that would cross a block boundary or not fit are skipped.

        Parameters
        ----------
        layout : list[ScheduledScanBlock]
            Contiguous science layout.
        placements : list[tuple[Time, str, Quantity, str]]
            Calibrator insertions as (start_time, block_name, duration, label).

        Returns
        -------
        list[ScheduledScanBlock]
            Combined science + calibrator schedule in chronological order.
        """
        if not placements:
            return list(layout)
        buckets: list[list[tuple[Time, str, u.Quantity, str]]] = [[] for _ in layout]
        for p in placements:
            ts = p[0]
            for bi, sb in enumerate(layout):
                if sb.start_time <= ts < sb.end_time:
                    buckets[bi].append(p)
                    break
        counters: dict[str, int] = {}
        result: list[ScheduledScanBlock] = []
        for bi, sb in enumerate(layout):
            plist = sorted(buckets[bi], key=lambda p: p[0].mjd)
            if not plist:
                result.append(sb)
                continue
            cur = sb.start_time
            for (ts, name, pdur, label) in plist:
                ts = max(ts, cur)
                te = ts + pdur
                if te > sb.end_time:
                    continue  # does not fit in the remaining block span
                if ts > cur:
                    chunk = self._make_sci_chunk(sb.name, cur, ts)
                    if chunk:
                        result.append(chunk)
                counters[label] = counters.get(label, 0) + 1
                result.append(self._make_aux_block(name, ts, te, f"{label}_{counters[label]}"))
                cur = te
            if cur < sb.end_time:
                chunk = self._make_sci_chunk(sb.name, cur, sb.end_time)
                if chunk:
                    result.append(chunk)
        return result

    def _schedule_cpsat(self, sci_names: list[str], ff_names: list[str],
                        t0: Time, t1: Time) -> list[ScheduledScanBlock]:
        """Constraint-programming schedule: science layout + CP-SAT FF/calibrator insertion.

        Builds a contiguous science layout, then inserts fringe finders (placed by
        CP-SAT), eMERLIN 3C286, and polcal scans by splitting the science blocks.

        Parameters
        ----------
        sci_names : list[str]
            Science block names.
        ff_names : list[str]
            FF block names.
        t0 : Time
            Observation start time.
        t1 : Time
            Observation end time.

        Returns
        -------
        list[ScheduledScanBlock]
            Full schedule in chronological order, or empty on failure.
        """
        eff_t0, eff_t1 = self._effective_obs_window(t0, t1)
        layout = self._layout_science(sci_names, eff_t0, eff_t1)
        if not layout:
            return []

        aux: list[tuple[Time, str, u.Quantity, str]] = []
        for ts, fn in self._cpsat_ff_placement(ff_names, layout, eff_t0, eff_t1):
            aux.append((ts, fn, self.FF_DUR, 'FF'))
        if self._has_emerlin():
            em = self._create_emerlin_3c286_block()
            if em:
                r = self._find_best(em, eff_t0, eff_t1, self.EMERLIN_3C286_DUR)
                if r:
                    aux.append((r[0], em, self.EMERLIN_3C286_DUR, 'eMERLIN_3C286'))
        if self._polcal:
            for ts, pc in self._polcal_target_times(eff_t0, eff_t1):
                aux.append((ts, pc, self.POLCAL_DUR, 'POLCAL'))

        return sorted(self._insert_aux(layout, aux), key=lambda b: b.start_time.mjd)

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def schedule(self) -> dict[str, ScanBlock]:
        """Generate the complete observation schedule.

        For a single science block: FF (evenly spread) → eMERLIN 3C286 → polcal → science fills gaps.
        For multiple science blocks: multi-source mode with equalised time per source,
        optimal window selection, slew minimisation, and FF scans distributed between
        source blocks.

        Returns
        -------
        dict[str, ScanBlock]
            Ordered mapping ``'NNN_label' → ScanBlock``.
        """
        t0, t1 = self.obs.times[0], self.obs.times[-1]

        sci_names = [n for n in self._sci_blocks()
                     if not n.startswith('FF_') and not n.startswith('POLCAL_') and n != 'eMERLIN_3C286']
        ff_names = self._select_fringefinders()

        # ---- Preferred path: CP-SAT constraint-programming layout ----
        if _HAS_ORTOOLS and sci_names:
            cp_sched = self._schedule_cpsat(sci_names, ff_names, t0, t1)
            if cp_sched:
                self._scheduled = sorted(cp_sched, key=lambda b: b.start_time.mjd)
                return {f"{i + 1:03d}_{b.name}": b.block for i, b in enumerate(self._scheduled)}
            log.warning("CP-SAT scheduling produced no blocks; falling back to heuristic scheduler.")

        if len(sci_names) > 1:
            # ---- Multi-source path ----
            all_sched: list[ScheduledScanBlock] = self._schedule_multi_source(sci_names, ff_names, t0, t1)

            # eMERLIN 3C286 and polcal are added in any remaining gaps
            if self._has_emerlin():
                em_block = self._create_emerlin_3c286_block()
                if em_block:
                    slots = self._get_slots(all_sched, t0, t1)
                    em = self._schedule_single_block(em_block, slots, 'eMERLIN_3C286')
                    if em:
                        all_sched.append(em)
            if self._polcal:
                polcal_names = self._create_polcal_blocks()
                slots = self._get_slots(all_sched, t0, t1)
                all_sched.extend(self._schedule_polcal(polcal_names, slots))
        else:
            # ---- Single-source path (original logic) ----
            ff_sched = self._schedule_ff(ff_names, t0, t1) if ff_names else []
            all_sched = list(ff_sched)

            if self._has_emerlin():
                em_block = self._create_emerlin_3c286_block()
                if em_block:
                    slots = self._get_slots(all_sched, t0, t1)
                    em = self._schedule_single_block(em_block, slots, 'eMERLIN_3C286')
                    if em:
                        all_sched.append(em)

            if self._polcal:
                polcal_names = self._create_polcal_blocks()
                slots = self._get_slots(all_sched, t0, t1)
                all_sched.extend(self._schedule_polcal(polcal_names, slots))

            slots = self._get_slots(all_sched, t0, t1)
            all_sched.extend(self._schedule_sci(sci_names, slots))

        self._scheduled = sorted(all_sched, key=lambda b: b.start_time.mjd)
        return {f"{i + 1:03d}_{b.name}": b.block for i, b in enumerate(self._scheduled)}

    def get_scheduled_blocks(self) -> list[ScheduledScanBlock]:
        """Return scheduled blocks in chronological order.

        Returns
        -------
        list[ScheduledScanBlock]
            Scheduled blocks sorted by start time.
        """
        return self._scheduled

    # ------------------------------------------------------------------
    # Key-file generation  (REWRITTEN — group/rep + Jb1)
    # ------------------------------------------------------------------

    @staticmethod
    def _scans_match(a: list[Scan], b: list[Scan]) -> bool:
        """Return True if two scan lists have the same sources and durations.

        Parameters
        ----------
        a : list[Scan]
        First scan list.
        b : list[Scan]
        Second scan list.

        Returns
        -------
        bool
        True if sources and durations match within 1 second tolerance.
        """
        if len(a) != len(b):
            return False
        return all(x.source.name == y.source.name
                   and abs(x.duration.to(u.s).value - y.duration.to(u.s).value) < 1.0
                   for x, y in zip(a, b))

    @staticmethod
    def _detect_cycle(scans: list[Scan]) -> tuple[list[Scan], int, list[Scan]]:
        """Find the repeating cycle in scans that covers the most total scans.

        Tries all candidate cycle lengths and picks the one where
        ``reps * cycle_len`` is maximised (i.e. covers the largest portion
        of the scan list).  This avoids the greedy trap of matching a short
        sub-cycle (e.g. pcal+target) when a longer cycle (including the
        check source every N) exists.

        Parameters
        ----------
        scans : list[Scan]
        List of scans to analyze.

        Returns
        -------
        tuple[list[Scan], int, list[Scan]]
        (cycle, n_reps, remainder) where cycle is one repetition of the pattern,
        n_reps is the number of full repetitions (≥ 1), and remainder is the
        left-over scans after the last full repetition.
        """
        n = len(scans)
        # Precomputed per-scan keys: integer source code and duration in seconds (same values as _scans_match).
        name_codes: dict[str, int] = {}
        codes = np.array([name_codes.setdefault(s.source.name, len(name_codes)) for s in scans], dtype=np.int64)
        durs = np.array([s.duration.to(u.s).value for s in scans], dtype=np.float64)
        best_clen, best_reps, best_covered = 0, 1, 0
        for clen in range(2, n // 2 + 1):
            n_chunks = n // clen
            m = n_chunks * clen
            ref = np.arange(m) % clen
            # Scan j matches the pattern when it equals scan (j mod clen), as in _scans_match.
            ok = (codes[:m] == codes[ref]) & (np.abs(durs[:m] - durs[ref]) < 1.0)
            chunk_ok = ok.reshape(n_chunks, clen).all(axis=1)
            reps = n_chunks if chunk_ok.all() else int(np.argmin(chunk_ok))
            covered = reps * clen
            if reps >= 2 and covered > best_covered:
                best_clen, best_reps, best_covered = clen, reps, covered
        if best_clen > 0:
            return scans[:best_clen], best_reps, scans[best_covered:]
        return scans, 1, []

    def _jb1_station_name(self) -> Optional[str]:
        """Return the SCHED-style station name for Jb1, or None.

        Returns
        -------
        str or None
        SCHED station name for JB1, or None if not in array.
        """
        for s in self.obs.stations:
            if s.codename.upper() == 'JB1':
                return s.name
        return None

    def _format_scan_line(self, scan: Scan, gap: str = '0:00', indent: str = '',
                          dur_override: Optional[u.Quantity] = None) -> str:
        """Format a single scan as a SCHED source= line.

        Parameters
        ----------
        scan : Scan
        Scan to format.
        gap : str, optional
        Gap before scan in 'M:SS' format. Default is '0:00'.
        indent : str, optional
        Indentation string. Default is ''.
        dur_override : Quantity or None, optional
        If given, use this duration instead of scan.duration.

        Returns
        -------
        str
        Formatted SCHED source line.
        """
        intent = _intent_str(scan.source.type)
        intent_part = f" intent='{intent}'" if intent else ''
        dur = dur_override if dur_override is not None else scan.duration
        return f"{indent}source='{_sched_safe(scan.source.name)}' gap={gap} dur={_fmt_dur(dur)}{intent_part} /"

    def _stations_line(self, exclude: Optional[set[str]] = None, indent: str = '') -> str:
        """Build a ``stations = ...`` line, optionally excluding codenames.

        Parameters
        ----------
        exclude : set[str] or None, optional
        Set of codenames to exclude. Default is None.
        indent : str, optional
        Indentation string. Default is ''.

        Returns
        -------
        str
        Formatted stations line using SCHED catalog names.
        """
        names = [s.sched_name for s in self.obs.stations
                 if exclude is None or s.codename.upper() not in exclude]
        return f"{indent}stations = {', '.join(names)}"

    def _build_science_scans(self, sb: ScheduledScanBlock) -> list[str]:
        """Convert a science ScheduledScanBlock into SCHED key-file lines.

        Detects repeating cycles and emits ``group N rep R``.
        If Jb1 is in the array, every other PHASECAL scan inside a group
        excludes Jb1 to stay under 12 source-changes / hour.

        Parameters
        ----------
        sb : ScheduledScanBlock
        Scheduled science block to convert.

        Returns
        -------
        list[str]
        List of SCHED key-file lines.
        """
        cycle, n_reps, remainder = self._detect_cycle(sb.scans)
        has_jb1 = self._has_jb1()
        all_stations = self._stations_line()
        lines: list[str] = []

        if n_reps >= 2:
            # Build one cycle with Jb1 handling
            cycle_lines = self._format_cycle(cycle, has_jb1, indent='    ')
            n_scan_lines = sum(1 for ln in cycle_lines if ln.strip().startswith("source="))
            lines.append(f"\ngroup {n_scan_lines} rep {n_reps}")
            lines.append(f"    {all_stations.strip()}")
            lines.extend(cycle_lines)
            # Remainder (partial last cycle) — also needs gaps
            if remainder:
                lines.append(all_stations)
                rem_gaps = self._compute_gap_positions(remainder)
                for idx, scan in enumerate(remainder):
                    gap_str = '0:00'
                    dur_ov: Optional[u.Quantity] = None
                    if idx in rem_gaps:
                        gap_str = _fmt_dur(self.GAP_DUR)
                        dur_ov = scan.duration - self.GAP_DUR
                    lines.append(self._format_scan_line(scan, gap=gap_str, dur_override=dur_ov))
        else:
            # No repeating cycle detected — still insert gaps
            all_gaps = self._compute_gap_positions(sb.scans)
            for idx, scan in enumerate(sb.scans):
                gap_str = '0:00'
                dur_ov = None
                if idx in all_gaps:
                    gap_str = _fmt_dur(self.GAP_DUR)
                    dur_ov = scan.duration - self.GAP_DUR
                lines.append(self._format_scan_line(scan, gap=gap_str, dur_override=dur_ov))
        return lines

    def _compute_gap_positions(self, cycle: list[Scan]) -> set[int]:
        """Determine which scan indices in the cycle should have a gap inserted.

        Distributes gaps evenly across the cycle at approximately GAP_INTERVAL
        spacing.  Only PHASECAL scans whose duration exceeds GAP_DUR are
        eligible.  The gap duration is later subtracted from the scan duration
        so the total time per scan slot stays the same.

        Parameters
        ----------
        cycle : list[Scan]
        List of scans in the cycle.

        Returns
        -------
        set[int]
        Indices into cycle that should get ``gap=0:30``.
        """
        total_dur = sum(s.duration.to(u.min).value for s in cycle)
        interval_min = self.GAP_INTERVAL.to(u.min).value
        gap_dur_s = self.GAP_DUR.to(u.s).value

        if total_dur < interval_min * 0.8:
            return set()

        # Find eligible phasecal positions with their cumulative start times
        pcal_positions: list[tuple[int, float]] = []
        cum = 0.0
        for i, scan in enumerate(cycle):
            if (scan.source.type == SourceType.PHASECAL
                    and scan.duration.to(u.s).value > gap_dur_s):
                pcal_positions.append((i, cum))
            cum += scan.duration.to(u.min).value

        if not pcal_positions:
            return set()

        # Compute how many gaps we need and their ideal positions
        n_gaps = max(1, round(total_dur / interval_min))
        n_gaps = min(n_gaps, len(pcal_positions))
        ideal_times = [(k + 0.5) * total_dur / n_gaps for k in range(n_gaps)]

        # Assign each ideal gap time to the nearest eligible phasecal
        used: set[int] = set()
        gap_positions: set[int] = set()
        for ideal_t in ideal_times:
            best_idx, best_dist = -1, float('inf')
            for scan_idx, cum_t in pcal_positions:
                if scan_idx in used:
                    continue
                dist = abs(cum_t - ideal_t)
                if dist < best_dist:
                    best_dist, best_idx = dist, scan_idx
            if best_idx >= 0:
                gap_positions.add(best_idx)
                used.add(best_idx)
        return gap_positions

    def _format_cycle(self, cycle: list[Scan], handle_jb1: bool, indent: str = '') -> list[str]:
        """Format one cycle of scans for the SCHED key file.

        - Inserts ~30 s gaps on phasecal scans every ~10 min of elapsed cycle
          time.  The gap is subtracted from the scan duration so the total time
          per scan slot stays the same.
        - Every other PHASECAL scan gets a station line excluding Jb1 when
          *handle_jb1* is True (to stay under 12 source-changes / hour).

        Parameters
        ----------
        cycle : list[Scan]
        List of scans in the cycle.
        handle_jb1 : bool
        Whether to handle Jb1 source-change mitigation.
        indent : str, optional
        Indentation string. Default is ''.

        Returns
        -------
        list[str]
        Formatted SCHED key-file lines.
        """
        lines: list[str] = []
        pcal_idx = 0
        all_stations_line = self._stations_line(indent=indent)
        no_jb1_line = self._stations_line(exclude={'JB1'}, indent=indent)
        gap_set = self._compute_gap_positions(cycle)

        for i, scan in enumerate(cycle):
            gap_str = '0:00'
            dur_override: Optional[u.Quantity] = None
            if i in gap_set:
                gap_str = _fmt_dur(self.GAP_DUR)
                dur_override = scan.duration - self.GAP_DUR

            if handle_jb1 and scan.source.type == SourceType.PHASECAL:
                pcal_idx += 1
                if pcal_idx % 2 == 0:
                    lines.append(no_jb1_line)
                    lines.append(self._format_scan_line(scan, gap=gap_str,
                                                        dur_override=dur_override, indent=indent))
                    lines.append(all_stations_line)
                    continue
            lines.append(self._format_scan_line(scan, gap=gap_str,
                                                dur_override=dur_override, indent=indent))
        return lines

    def generate_key_file(self, experiment_code: str = 'EXCODE', pi_name: str = 'PI Name',
                          pi_email: str = 'pi@example.com', pi_institute: str = 'Institute',
                          setup_file: Optional[str] = None, comments: str = '',
                          template_path: Optional[str] = None) -> str:
        """Generate a SCHED .key file from the current schedule.

        Uses ``group N rep R`` for repeated science cycles and excludes Jb1
        from every other phase-cal scan when Jb1 is present. User-provided strings
        (experiment code, PI fields, comments, source names) have newlines and single
        quotes stripped before insertion so they cannot break out of SCHED fields.

        Parameters
        ----------
        experiment_code : str
            Experiment code for the observation.
        pi_name : str
            Principal investigator name.
        pi_email : str
            Principal investigator email.
        pi_institute : str
            Principal investigator institute.
        setup_file : str or None, optional
            Frequency setup file name. If None, a setup file is guessed from the
            array, band and data rate by searching ``~/.pysched/setups``.
        comments : str, optional
            Additional comments for the key file. Default is ''.
        template_path : str or None, optional
            Path to a custom SCHED .key template. If None, the bundled template
            is used.

        Returns
        -------
        str
            Complete .key file content.

        Raises
        ------
        ValueError
            If the template is missing the mandatory ``{SCANS}``, ``{SOURCES}``,
            or start-of-observation placeholders, or if their generated values
            are empty.
        """
        if template_path is None:
            template_file = resources.files('vlbiplanobs.data').joinpath('key_file.key.template')
        else:
            template_file = Path(template_path)
        template = Path(template_file).read_text(encoding='utf-8')

        mandatory_placeholders = ('{SCANS}', '{SOURCES}', '{YEAR}', '{MONTH}', '{DAY}', '{START_TIME}')
        missing = [p for p in mandatory_placeholders if p not in template]
        if missing:
            raise ValueError(f"Template is missing mandatory placeholder(s): {', '.join(missing)}")

        # ---- Source catalog ----
        all_sources: dict[str, Source] = {}
        for sb in self._scheduled:
            for scan in sb.scans:
                if scan.source.name not in all_sources:
                    all_sources[scan.source.name] = scan.source
        src_lines = []
        for src in all_sources.values():
            ra = src.coord.ra.to_string(unit=u.hourangle, sep=':', precision=4, pad=True)
            dec = src.coord.dec.to_string(unit=u.degree, sep=':', precision=3, pad=True, alwayssign=True)
            src_lines.append(f"  source='{_sched_safe(src.name)}' ra={ra} dec={dec} equinox='J2000' /")

        if not src_lines:
            raise ValueError("No sources were scheduled; {SOURCES} cannot be empty.")

        # ---- Scan section ----
        all_stations = self._stations_line()
        scan_lines: list[str] = [all_stations, '']
        for sb in self._scheduled:
            is_ff = sb.block.has(SourceType.FRINGEFINDER)
            is_emerlin = sb.block.has(SourceType.AMPLITUDECAL)
            is_polcal = sb.block.has(SourceType.POLCAL)

            if is_ff or is_emerlin or is_polcal:
                gap = '1:30' if is_emerlin else '0:00'
                for scan in sb.scans:
                    scan_lines.append(self._format_scan_line(scan, gap=gap))
            else:
                scan_lines.extend(self._build_science_scans(sb))
            scan_lines.append('')  # blank line between blocks

        if not self._scheduled:
            raise ValueError("No scans were scheduled; {SCANS} cannot be empty.")

        # ---- Setup ----
        if setup_file is None:
            setup_file = guess_setup_file(self.obs)
        setup_str = format_setup_line(setup_file)

        obs_mode = (f"{self.obs.band} {int(self.obs.datarate.to(u.Mbit / u.s).value)} Mbps"
                    if self.obs.band and self.obs.datarate is not None else "VLBI")
        start_time = self.obs.times[0]
        replacements = {
            'GENERATION_DATE': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'EXPERIMENT_CODE': _sched_safe(experiment_code).upper(),
            'PI_NAME': _sched_safe(pi_name), 'PI_EMAIL': _sched_safe(pi_email),
            'PI_INSTITUTE': _sched_safe(pi_institute), 'OBS_MODE': obs_mode, 'COMMENTS': _sched_safe(comments),
            'CORAVG': str(int(self.obs.inttime.to(u.s).value)),
            'CORCHAN': str(self.obs.channels) if self.obs.channels else '32',
            'CORNANT': str(len(self.obs.stations)),
            'STATIONS_CATALOG': 'none',
            'SOURCES': '\n'.join(src_lines),
            'SETUP': setup_str,
            'YEAR': str(start_time.datetime.year),
            'MONTH': str(start_time.datetime.month),
            'DAY': str(start_time.datetime.day),
            'START_TIME': start_time.datetime.strftime('%H:%M:%S'),
            'STATIONS': ', '.join(s.sched_name for s in self.obs.stations),
            'SCANS': '\n'.join(scan_lines),
        }

        def _replace_placeholder(match: re.Match) -> str:
            key = match.group(1)
            return str(replacements.get(key, match.group(0)))

        result = re.sub(r'\{([A-Za-z_][A-Za-z0-9_]*)\}', _replace_placeholder, template)
        return result
