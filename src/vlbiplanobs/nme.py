# -*- coding: utf-8 -*-
# Licensed under GPLv3+ - see LICENSE
"""Network Monitoring Experiment (NME) planning and SCHED key-file generation.

An NME is a test observation where the full time is covered by fringe-finder
scans (~15 min each) that all participating antennas can observe. At regular
intervals, a scan carries a SCHED ``grabto='FILE' grabtime=...`` instruction so
the stations transfer a small chunk of data (ftp fringe test) to the correlator.

Grab cadence (all times relative to the observation start):
- duration > 2.5 h: first grab at +10 min, then every 30 min.
- duration <= 2.5 h: first grab at +5 min, then every 15 min, plus one grab
  right before the end of the observation.

SCHED semantics (GRABTIME): the first value is the number of seconds of data to
send; the second is how many seconds before the end of the scan the grabbed data
end. With ``grabtime=2,118`` the grabbed data start 120 s before the scan end, so
the grab scan ends GRAB_DATA_S + GRAB_END_OFFSET_S after the nominal grab time.
"""

from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Optional
import logging
import re
import numpy as np
from astropy import units as u
from astropy.time import Time
from .sources import Source, SourceType
from .stations import Stations
from . import calibrators

log = logging.getLogger(__name__)

GRAB_DATA_S = 2
GRAB_END_OFFSET_S = 118
GRAB_TAIL_S = GRAB_DATA_S + GRAB_END_OFFSET_S
LONG_OBS_THRESHOLD_S = int(2.5 * 3600)
LONG_FIRST_GRAB_S, LONG_GRAB_STEP_S = 10 * 60, 30 * 60
SHORT_FIRST_GRAB_S, SHORT_GRAB_STEP_S = 5 * 60, 15 * 60
MIN_GRAB_SEPARATION_S = 5 * 60
MIN_DURATION_S = 15 * 60
TARGET_SCAN_S = 15 * 60
SAMPLE_STEP_S = 60
DEFAULT_MIN_ELEVATION = 15.0 * u.deg
DEFAULT_MIN_FLUX = 1.0 * u.Jy


@dataclass
class NMEScan:
    """One scan of the NME. Offsets are integer seconds from the observation start.

    ``slot_start_s`` is when the previous scan ends; recording starts at
    ``slot_start_s + gap_s`` and ends at ``stop_s``.
    """
    slot_start_s: int
    gap_s: int
    stop_s: int
    grab_s: Optional[int] = None
    source: Optional[Source] = None
    n_visible: int = 0

    @property
    def rec_start_s(self) -> int:
        """Second (from observation start) when on-source recording starts."""
        return self.slot_start_s + self.gap_s

    @property
    def dur_s(self) -> int:
        """Recording duration in seconds."""
        return self.stop_s - self.rec_start_s


def grab_offsets(duration_s: int) -> list[int]:
    """Return the nominal grab times (seconds from start) for an NME of the given duration.

    Parameters
    ----------
    duration_s : int
        Total duration of the observation in seconds.

    Returns
    -------
    list[int]
        Sorted grab times. Each grab scan ends GRAB_TAIL_S after its grab time.

    Raises
    ------
    ValueError
        If the observation is shorter than MIN_DURATION_S.
    """
    if duration_s < MIN_DURATION_S:
        raise ValueError(f"An NME must last at least {MIN_DURATION_S // 60} min (got {duration_s / 60:.1f} min).")

    if duration_s > LONG_OBS_THRESHOLD_S:
        return list(range(LONG_FIRST_GRAB_S, duration_s - GRAB_TAIL_S + 1, LONG_GRAB_STEP_S))

    final = duration_s - GRAB_TAIL_S
    regular = [g for g in range(SHORT_FIRST_GRAB_S, final + 1, SHORT_GRAB_STEP_S)
               if g <= final - MIN_GRAB_SEPARATION_S]
    return regular + [final]


def _split_segment(length_s: int, n_scans: int) -> list[int]:
    """Split a segment into n slot lengths of whole minutes; the last slot takes the remainder."""
    base = (length_s // n_scans) // 60 * 60
    return [base] * (n_scans - 1) + [length_s - base * (n_scans - 1)]


def _gap_for_slot(slot_s: int) -> int:
    """Return the slewing gap (s) placed at the beginning of a scan slot of the given length."""
    if slot_s >= 12 * 60:
        return 5 * 60
    if slot_s >= 5 * 60:
        return 2 * 60
    return 0


def plan_scans(duration_s: int, grabs: list[int]) -> list[NMEScan]:
    """Cover the full observation with ~15-min scans, aligning grab scans to the grab times.

    Each grab scan ends exactly GRAB_TAIL_S after its grab time. The time between
    consecutive grab-scan ends is split into round(L / 15 min) scans (at least one;
    the first segment gets at least two when >= 10 min, so the first scan is not a
    grab scan). The first scan of the observation has no gap.

    Parameters
    ----------
    duration_s : int
        Total duration of the observation in seconds.
    grabs : list[int]
        Grab times in seconds from the start (see ``grab_offsets``).

    Returns
    -------
    list[NMEScan]
        Contiguous scans covering [0, duration_s], without sources assigned yet.
    """
    grab_by_end = {g + GRAB_TAIL_S: g for g in grabs}
    boundaries = sorted({0, duration_s, *grab_by_end.keys()})
    scans: list[NMEScan] = []
    for seg_idx, (a, b) in enumerate(zip(boundaries[:-1], boundaries[1:])):
        length = b - a
        n_scans = max(1, round(length / TARGET_SCAN_S))
        if seg_idx == 0 and length >= 10 * 60:
            n_scans = max(n_scans, 2)
        t = a
        for slot in _split_segment(length, n_scans):
            gap = 0 if not scans else _gap_for_slot(slot)
            scans.append(NMEScan(slot_start_s=t, gap_s=gap, stop_s=t + slot))
            t += slot
        if b in grab_by_end:
            scans[-1].grab_s = grab_by_end[b]
            if scans[-1].dur_s < GRAB_TAIL_S + 60:
                scans[-1].gap_s = 0
    return scans


def resolve_fringe_finders(names: list[str]) -> list[Source]:
    """Resolve user-given fringe-finder names into sources (RFC catalog first, then name/coords lookup).

    Parameters
    ----------
    names : list[str]
        Source names, or 'name/coordinates' strings.

    Returns
    -------
    list[Source]
        Resolved sources.

    Raises
    ------
    ValueError
        If a source cannot be resolved.
    """
    catalog = calibrators.RFCCatalog(min_flux=0.0 * u.Jy, band='c', include_missing=True)
    resolved: list[Source] = []
    for name in names:
        src = catalog.get_source(name) if '/' not in name else None
        if src is None:
            try:
                src = Source.source_from_str(name, source_type=SourceType.FRINGEFINDER)
            except Exception as error:
                raise ValueError(f"Fringe finder '{name}' could not be resolved: {error}")
        resolved.append(src)
    return resolved


def candidate_fringe_finders(band: str, min_flux: u.Quantity = DEFAULT_MIN_FLUX) -> list[Source]:
    """Return RFC sources with an unresolved flux above ``min_flux`` at the given band.

    Parameters
    ----------
    band : str
        Observing band (e.g. '18cm').
    min_flux : Quantity
        Minimum unresolved flux density.

    Returns
    -------
    list[Source]
        Candidate fringe finders sorted by decreasing unresolved flux at the band.
    """
    rfc_band = calibrators._wavelength_to_rfc_band(band)
    catalog = calibrators.RFCCatalog(min_flux=min_flux, band=rfc_band)
    sources = list(catalog.sources)
    sources.sort(key=lambda s: -source_flux(s, band))
    log.info("NME: %d RFC fringe-finder candidates above %s at band %s", len(sources), min_flux, band)
    return sources


def source_flux(src: Source, band: str) -> float:
    """Unresolved flux (Jy) of a source at the band, or 0.0 for non-RFC sources."""
    if isinstance(src, calibrators.CalibratorSource):
        return src.get_flux_at_band(band)[1]
    return 0.0


def visibility_counts(stations: Stations, sources: list[Source], times: Time,
                      min_elevation: u.Quantity = DEFAULT_MIN_ELEVATION) -> np.ndarray:
    """Count how many stations can observe each source at each time.

    Parameters
    ----------
    stations : Stations
        Participating stations.
    sources : list[Source]
        Sources to evaluate.
    times : Time
        Sampling times.
    min_elevation : Quantity
        Minimum elevation required (on top of each station's mount limits).

    Returns
    -------
    np.ndarray
        Integer array of shape (n_times, n_sources).
    """
    ra_rad = np.array([s.coord.ra.rad for s in sources], dtype=np.float64)
    dec_rad = np.array([s.coord.dec.rad for s in sources], dtype=np.float64)
    dec_deg = np.degrees(dec_rad)
    min_el = min_elevation.to(u.deg).value
    counts = np.zeros((len(times), len(sources)), dtype=np.int32)
    for station in stations:
        elev, az, ha_hours = calibrators._batch_altaz_erfa(ra_rad, dec_rad, times, station)
        mask = calibrators._station_observable_mask(elev, az, ha_hours, dec_deg, station) & (elev >= min_el)
        counts += mask.astype(np.int32)
    return counts


def _scan_min_counts(scan: NMEScan, counts: np.ndarray) -> np.ndarray:
    """Minimum number of visible stations per source over the scan recording interval."""
    i0 = scan.rec_start_s // SAMPLE_STEP_S
    i1 = -(-scan.stop_s // SAMPLE_STEP_S)
    return counts[i0:i1 + 1].min(axis=0)


def assign_sources(scans: list[NMEScan], sources: list[Source], counts: np.ndarray, n_stations: int) -> None:
    """Assign a fringe finder to each scan, in place.

    Prefers sources observable by all stations during the whole scan, keeps the
    previous source while possible (fewer slews), otherwise picks the source that
    stays fully visible for the most consecutive scans (ties: order in ``sources``,
    i.e. brightest first). If no source is visible by all stations, the one seen
    by the most stations is used and a warning is logged.

    Parameters
    ----------
    scans : list[NMEScan]
        Scans to fill (``source`` and ``n_visible`` are set).
    sources : list[Source]
        Candidate sources, in order of preference.
    counts : np.ndarray
        Output of ``visibility_counts`` sampled every SAMPLE_STEP_S from the start.
    n_stations : int
        Number of participating stations.
    """
    scan_counts = np.array([_scan_min_counts(s, counts) for s in scans])  # (n_scans, n_sources)
    full = scan_counts >= n_stations
    current: Optional[int] = None
    for k, scan in enumerate(scans):
        if current is not None and full[k, current]:
            choice = current
        elif full[k].any():
            run = np.zeros(len(sources), dtype=np.int32)
            alive = full[k].copy()
            for j in range(k, len(scans)):
                alive &= full[j]
                if not alive.any():
                    break
                run += alive
            choice = int(np.argmax(run))
        else:
            choice = int(np.argmax(scan_counts[k]))
            log.warning("NME scan %d: no source visible by all %d stations; using %s (%d stations)",
                        k + 1, n_stations, sources[choice].name, scan_counts[k, choice])
        scan.source = sources[choice]
        scan.n_visible = int(scan_counts[k, choice])
        current = choice


def plan_nme(stations: Stations, start: Time, duration: u.Quantity, band: str,
             fringefinders: Optional[list[Source]] = None,
             min_elevation: u.Quantity = DEFAULT_MIN_ELEVATION,
             min_flux: u.Quantity = DEFAULT_MIN_FLUX) -> tuple[list[NMEScan], list[Source], np.ndarray, Time]:
    """Plan a full NME: grab times, scans, and fringe-finder assignment.

    Parameters
    ----------
    stations : Stations
        Participating stations.
    start : Time
        Start of the observation (UTC).
    duration : Quantity
        Duration of the observation.
    band : str
        Observing band (e.g. '18cm').
    fringefinders : list[Source] or None
        Fringe finders to use. If None, RFC sources brighter than ``min_flux`` are considered.
    min_elevation : Quantity
        Minimum elevation for a station to count as observing a source.
    min_flux : Quantity
        Minimum unresolved flux for automatic candidates.

    Returns
    -------
    tuple[list[NMEScan], list[Source], np.ndarray, Time]
        (scans, candidate sources, visibility counts (n_times, n_sources), sample times).

    Raises
    ------
    ValueError
        If the duration is too short or no candidate fringe finders exist.
    """
    duration_s = int(round(duration.to(u.s).value))
    grabs = grab_offsets(duration_s)
    scans = plan_scans(duration_s, grabs)
    sources = fringefinders if fringefinders else candidate_fringe_finders(band, min_flux)
    if not sources:
        raise ValueError(f"No fringe-finder candidates found above {min_flux} at {band}.")
    times = start + np.arange(0, duration_s + SAMPLE_STEP_S, SAMPLE_STEP_S) * u.s
    counts = visibility_counts(stations, sources, times, min_elevation)
    assign_sources(scans, sources, counts, len(stations))
    log.info("NME planned: %d scans, %d grabs, sources used: %s", len(scans), len(grabs),
             sorted({s.source.name for s in scans}))
    return scans, sources, counts, times


def visible_windows(mask: np.ndarray, times: Time) -> list[tuple[Time, Time]]:
    """Convert a boolean time mask into a list of contiguous (start, end) time windows."""
    windows: list[tuple[Time, Time]] = []
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return windows
    breaks = np.flatnonzero(np.diff(idx) > 1)
    starts = np.concatenate(([idx[0]], idx[breaks + 1]))
    ends = np.concatenate((idx[breaks], [idx[-1]]))
    return [(times[a], times[b]) for a, b in zip(starts, ends)]


def _fmt_ms(seconds: int) -> str:
    """Format seconds as SCHED 'M:SS'."""
    return f"{seconds // 60}:{seconds % 60:02d}"


def _source_catalog_line(src: Source) -> str:
    """Build a SCHED srccat line for a source (J2000 name plus IVS alias when known)."""
    ra = src.coord.ra.to_string(unit=u.hourangle, sep=':', precision=6, pad=True)
    dec = src.coord.dec.to_string(unit=u.degree, sep=':', precision=5, pad=True, alwayssign=True)
    ivs = getattr(src, 'ivsname', None)
    names = f"'{src.name}', '{ivs}'" if ivs and ivs != src.name else f"'{src.name}'"
    return f"source={names} ra={ra} dec={dec} equinox='J2000' /"


def scan_lines(scans: list[NMEScan], start: Time, n_stations: int) -> list[str]:
    """Render the NME scans as SCHED key-file lines (grab scans toggle grabto FILE/NONE)."""
    lines: list[str] = []
    grab_active = False
    for k, scan in enumerate(scans, start=1):
        if scan.n_visible < n_stations:
            lines.append(f"! WARNING: only {scan.n_visible}/{n_stations} antennas can observe this scan.")
        if scan.grab_s is not None:
            grab_utc = (start + scan.grab_s * u.s).datetime.strftime('%H:%M:%S')
            lines.append(f"! grab: scan {k} at {grab_utc}")
            lines.append(f"grabto='FILE' grabtime={GRAB_DATA_S},{GRAB_END_OFFSET_S}")
            grab_active = True
        elif grab_active:
            lines.append("grabto='NONE'")
            grab_active = False
        lines.append(f"source='{scan.source.name}' gap={_fmt_ms(scan.gap_s)} dur={_fmt_ms(scan.dur_s)} /")
        lines.append('')
    return lines


def grab_summary(scans: list[NMEScan], start: Time) -> str:
    """Cover-letter lines listing each ftp fringe test: 'HH:MM:SS (scan  N, 2 sec, SOURCE)'."""
    rows = [f"{(start + s.grab_s * u.s).datetime.strftime('%H:%M:%S')} (scan {k:2d}, {GRAB_DATA_S} sec, "
            f"{s.source.name})" for k, s in enumerate(scans, start=1) if s.grab_s is not None]
    return '\n'.join(rows)


def generate_nme_key_file(scans: list[NMEScan], stations: Stations, start: Time, band: str,
                          experiment_code: str, setup_line: str, datarate_mbps: Optional[int] = None,
                          pi_name: str = 'PI Name', pi_email: str = 'pi@example.com',
                          pi_institute: str = 'Institute', phone: str = '',
                          template_path: Optional[str] = None) -> str:
    """Generate the SCHED .key file content for an NME.

    Parameters
    ----------
    scans : list[NMEScan]
        Planned scans with sources assigned (see ``plan_nme``).
    stations : Stations
        Participating stations.
    start : Time
        Start of the observation (UTC).
    band : str
        Observing band (e.g. '18cm').
    experiment_code : str
        Experiment code (e.g. 'N25L1').
    setup_line : str
        Value for the ``setup = ...`` line (see ``scheduler.format_setup_line``).
    datarate_mbps : int or None
        Data rate written in the notes.
    pi_name, pi_email, pi_institute, phone : str
        Cover information.
    template_path : str or None
        Custom template. If None, the bundled ``nme_key_file.key.template`` is used.

    Returns
    -------
    str
        Complete .key file content. Unknown template placeholders are left untouched.

    Raises
    ------
    ValueError
        If the template lacks a mandatory placeholder or no scans are given.
    """
    if not scans or any(s.source is None for s in scans):
        raise ValueError("All NME scans must have a source assigned before writing the key file.")
    if template_path is None:
        template = resources.files('vlbiplanobs.data').joinpath('nme_key_file.key.template').read_text('utf-8')
    else:
        template = Path(template_path).read_text(encoding='utf-8')
    mandatory = ('{SCANS}', '{SOURCES}', '{YEAR}', '{MONTH}', '{DAY}', '{START_TIME}')
    missing = [p for p in mandatory if p not in template]
    if missing:
        raise ValueError(f"Template is missing mandatory placeholder(s): {', '.join(missing)}")

    used: dict[str, Source] = {}
    for scan in scans:
        used.setdefault(scan.source.name, scan.source)
    n_stations = len(stations)
    start_dt = start.datetime
    replacements = {
        'EXPERIMENT_CODE': experiment_code.upper(), 'BAND_LABEL': band,
        'DATARATE_MBPS': str(datarate_mbps) if datarate_mbps is not None else 'N/A',
        'PI_NAME': pi_name, 'PI_EMAIL': pi_email, 'PI_INSTITUTE': pi_institute, 'PHONE': phone,
        'DATE_LONG': f"{start_dt.day} {start_dt.strftime('%B %Y')}",
        'N_STATIONS': str(n_stations), 'CORNANT': str(n_stations),
        'STATION_CODES': ', '.join(s.codename for s in stations),
        'GRAB_SUMMARY': grab_summary(scans, start),
        'SOURCES': '\n'.join(_source_catalog_line(s) for s in used.values()),
        'SETUP': setup_line,
        'YEAR': str(start_dt.year), 'MONTH': f"{start_dt.month:02d}", 'DAY': f"{start_dt.day:02d}",
        'START_TIME': start_dt.strftime('%H:%M:%S'),
        'STATIONS': ', '.join(s.sched_name for s in stations),
        'SCANS': '\n'.join(scan_lines(scans, start, n_stations)).rstrip(),
    }
    return re.sub(r'\{([A-Za-z_][A-Za-z0-9_]*)\}',
                  lambda m: str(replacements.get(m.group(1), m.group(0))), template)
