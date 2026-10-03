"""Server-side validation of every user-controlled observation parameter in the GUI.

Browser-side limits (Input min/max, dropdown options, dcc.Upload max_size) are NOT enforced
on the server: any value that reaches a Dash callback (direct callback POSTs, the shared-link
``?config=`` payload parsed by ``url_open``, localStorage-backed stores, ``store-obs-params``)
is attacker-controlled. Every value must pass through this module before it is used in a
computation.

All limits are declared as module-level constants so the UI (inputs.py) and the server checks
share one source of truth.
"""
from __future__ import annotations
import math
from datetime import datetime as dt
from typing import Any, Optional
from loguru import logger
from astropy.time import Time
from vlbiplanobs import freqsetups as fs
from vlbiplanobs import observation

# ----------------------------------------------------------------------------------------
# Limits (single source of truth; inputs.py reuses MAX_DURATION_H / MAX_UPLOAD_BYTES).
# ----------------------------------------------------------------------------------------
MAX_DURATION_H: float = 50.0       # Same as the 'duration' dbc.Input max in inputs.py.
MAX_TARGETS: int = 20              # Maximum number of target sources computed at once.
MAX_TARGET_LEN: int = 80           # Same as the limit in callbacks._validate_source_spec.
MAX_UPLOAD_BYTES: int = 50_000     # Maximum size of an uploaded source-list file.
MIN_ONSOURCE_PCT: float = 20.0     # Same as the 'onsourcetime' slider min.
MAX_ONSOURCE_PCT: float = 100.0    # Same as the 'onsourcetime' slider max.
MIN_START_YEAR: int = 1950         # Same as the 'startdate' DatePickerSingle min_date_allowed.
MAX_START_YEAR: int = 2100         # Same as the 'startdate' DatePickerSingle max_date_allowed.

DEFAULT_DATARATE: int = 2048
DEFAULT_SUBBANDS: int = 8
DEFAULT_CHANNELS: int = 64
DEFAULT_POLARIZATIONS: int = 2
DEFAULT_INTTIME: float = 2
DEFAULT_ONTARGET: float = 0.7


class InvalidObsParams(ValueError):
    """Raised when a user-provided observation parameter is out of the allowed values."""


def _is_number(value: Any) -> bool:
    """Return True if value is a finite int/float (bools are rejected)."""
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def band_name_from_index(band_index: Any) -> Optional[str]:
    """Convert the band-slider index into a band name, validating it.

    Parameters
    ----------
    band_index : Any
        Value of the 'band-slider' component. 0 or None means no band selected.

    Returns
    -------
    str or None
        Band name (a key of freqsetups.bands), or None when no band is selected.

    Raises
    ------
    InvalidObsParams
        If the index is not an integer within [0, len(freqsetups.bands)].
    """
    if band_index is None or band_index == 0:
        return None
    if not isinstance(band_index, int) or isinstance(band_index, bool) or not 0 < band_index <= len(fs.bands):
        raise InvalidObsParams(f"Invalid band index: {band_index!r}.")
    return list(fs.bands.keys())[band_index - 1]


def clean_target_specs(target_specs: Any) -> list[str]:
    """Return a safe list of target specs: strings only, stripped, non-empty, deduplicated.

    Entries longer than MAX_TARGET_LEN are dropped and the list is truncated to MAX_TARGETS.

    Parameters
    ----------
    target_specs : Any
        Raw value (expected list[str]) coming from store-targets or store-obs-params.

    Returns
    -------
    list[str]
        Cleaned list with at most MAX_TARGETS entries (order preserved).
    """
    if not isinstance(target_specs, (list, tuple)):
        return []
    cleaned: list[str] = []
    for spec in target_specs:
        if not isinstance(spec, str):
            continue
        spec = spec.strip()
        if not spec or len(spec) > MAX_TARGET_LEN or spec in cleaned:
            continue
        cleaned.append(spec)
        if len(cleaned) >= MAX_TARGETS:
            break
    if len(cleaned) >= MAX_TARGETS and len(target_specs) > MAX_TARGETS:
        logger.warning(f"Target list truncated to {MAX_TARGETS} entries (received {len(target_specs)}).")
    return cleaned


def clean_station_codenames(codenames: Any) -> list[str]:
    """Return the known station codenames from a raw list (unknown/non-string entries dropped).

    Parameters
    ----------
    codenames : Any
        Raw value (expected list[str]) of selected station codenames.

    Returns
    -------
    list[str]
        Sorted, deduplicated list of codenames present in observation._STATIONS.
    """
    if not isinstance(codenames, (list, tuple)):
        return []
    return sorted({c for c in codenames if isinstance(c, str) and c in observation._STATIONS})


def validate_duration(duration: Any) -> Optional[float]:
    """Validate the observing duration in hours.

    Parameters
    ----------
    duration : Any
        Raw duration value. None or <= 0 means "no duration given".

    Returns
    -------
    float or None
        Duration in hours, or None when not given.

    Raises
    ------
    InvalidObsParams
        If not a finite number or larger than MAX_DURATION_H.
    """
    if duration is None:
        return None
    if not _is_number(duration):
        raise InvalidObsParams(f"Duration must be a number of hours (got {duration!r}).")
    if duration <= 0:
        return None
    if duration > MAX_DURATION_H:
        raise InvalidObsParams(f"Duration must be at most {MAX_DURATION_H:g} hours (got {duration:g}).")
    return float(duration)


def ontarget_from_percent(percent: Any) -> float:
    """Convert the on-source time percentage into a fraction, validating it.

    Parameters
    ----------
    percent : Any
        Raw 'onsourcetime' slider value (percentage). None or 0 gives DEFAULT_ONTARGET.

    Returns
    -------
    float
        On-target fraction in (0, 1].

    Raises
    ------
    InvalidObsParams
        If not a number within [MIN_ONSOURCE_PCT, MAX_ONSOURCE_PCT].
    """
    if percent is None or percent == 0:
        return DEFAULT_ONTARGET
    if not _is_number(percent) or not MIN_ONSOURCE_PCT <= percent <= MAX_ONSOURCE_PCT:
        raise InvalidObsParams(f"On-source time must be between {MIN_ONSOURCE_PCT:g} and "
                               f"{MAX_ONSOURCE_PCT:g} % (got {percent!r}).")
    return percent / 100


def _validate_ontarget_fraction(ontarget: Any) -> float:
    """Validate an on-target fraction (as stored in store-obs-params).

    Raises
    ------
    InvalidObsParams
        If not a number within [MIN_ONSOURCE_PCT/100, 1].
    """
    if ontarget is None:
        return DEFAULT_ONTARGET
    if not _is_number(ontarget) or not MIN_ONSOURCE_PCT / 100 <= ontarget <= MAX_ONSOURCE_PCT / 100:
        raise InvalidObsParams(f"Invalid on-target fraction: {ontarget!r}.")
    return float(ontarget)


def _validate_choice(name: str, value: Any, allowed: dict, default: Any) -> Any:
    """Validate that value is one of the keys of a freqsetups dictionary.

    Parameters
    ----------
    name : str
        Parameter name used in the error message.
    value : Any
        Raw value. None gives default. Numeric strings are accepted (Dropdown values may be str).
    allowed : dict
        freqsetups dictionary whose keys are the allowed values.
    default : Any
        Value returned when value is None.

    Returns
    -------
    Any
        The matching allowed key.

    Raises
    ------
    InvalidObsParams
        If value is not one of the allowed keys.
    """
    if value is None:
        return default
    if isinstance(value, str):
        try:
            value = float(value)
        except ValueError:
            raise InvalidObsParams(f"Invalid {name}: {value!r}.")
    if not _is_number(value):
        raise InvalidObsParams(f"Invalid {name}: {value!r}.")
    for key in allowed:
        if key == value:
            return key
    raise InvalidObsParams(f"Invalid {name}: {value!r}. Allowed: {', '.join(str(k) for k in allowed)}.")


def validate_datarate(value: Any) -> int:
    """Validate the data rate (Mbit/s) against freqsetups.data_rates. See _validate_choice."""
    return _validate_choice('data rate', value, fs.data_rates, DEFAULT_DATARATE)


def validate_subbands(value: Any) -> int:
    """Validate the number of subbands against freqsetups.subbands. See _validate_choice."""
    return _validate_choice('number of subbands', value, fs.subbands, DEFAULT_SUBBANDS)


def validate_channels(value: Any) -> int:
    """Validate the number of channels against freqsetups.channels. See _validate_choice."""
    return _validate_choice('number of channels', value, fs.channels, DEFAULT_CHANNELS)


def validate_polarizations(value: Any) -> int:
    """Validate the number of polarizations against freqsetups.polarizations. See _validate_choice."""
    return _validate_choice('number of polarizations', value, fs.polarizations, DEFAULT_POLARIZATIONS)


def validate_inttime(value: Any) -> float:
    """Validate the integration time (s) against freqsetups.inttimes. See _validate_choice."""
    return _validate_choice('integration time', value, fs.inttimes, DEFAULT_INTTIME)


def validate_band(band: Any) -> str:
    """Validate a band name.

    Raises
    ------
    InvalidObsParams
        If band is not a key of freqsetups.bands.
    """
    if not isinstance(band, str) or band not in fs.bands:
        raise InvalidObsParams(f"Unknown band: {band!r}.")
    return band


def parse_start_time(startdate: Any, starttime: Any) -> Time:
    """Parse the start date ('YYYY-MM-DD', optionally followed by a time part) and time ('HH:MM').

    Parameters
    ----------
    startdate : Any
        Raw 'startdate' DatePickerSingle value.
    starttime : Any
        Raw 'starttime' Dropdown value.

    Returns
    -------
    astropy.time.Time
        Start time in UTC.

    Raises
    ------
    InvalidObsParams
        If either value is not a string, cannot be parsed, or the year is out of range.
    """
    if not isinstance(startdate, str) or not isinstance(starttime, str):
        raise InvalidObsParams("Start date and time must be strings.")
    try:
        start = dt.strptime(f"{startdate[:10]} {starttime}", '%Y-%m-%d %H:%M')
    except ValueError:
        raise InvalidObsParams(f"Invalid start date/time: {startdate[:20]!r} {starttime[:10]!r}.")
    if not MIN_START_YEAR <= start.year <= MAX_START_YEAR:
        raise InvalidObsParams(f"Start year must be within {MIN_START_YEAR}-{MAX_START_YEAR}.")
    return Time(start, format='datetime', scale='utc')


def normalize_obs_params(raw: Any) -> dict:
    """Validate and normalize a serialized observation-parameters dict (store-obs-params shape).

    Expected keys: band (str), stations (list[str]), targets (list[str] or None),
    duration (h, or None), ontarget (fraction), startdate/starttime (str or None),
    datarate (Mbit/s), subbands, channels, polarizations, inttime (s).

    Parameters
    ----------
    raw : Any
        Raw dictionary.

    Returns
    -------
    dict
        Same keys with validated, normalized values; 'targets' is None when no valid targets and
        an extra 'start_time' key holds the parsed astropy Time (or None).

    Raises
    ------
    InvalidObsParams
        If raw is not a dict or any parameter is invalid, or there are no valid stations.
    """
    if not isinstance(raw, dict):
        raise InvalidObsParams("Observation parameters must be a dictionary.")
    stations = clean_station_codenames(raw.get('stations'))
    if not stations:
        raise InvalidObsParams("No valid stations selected.")
    targets = clean_target_specs(raw.get('targets'))
    startdate, starttime = raw.get('startdate'), raw.get('starttime')
    start_time = parse_start_time(startdate, starttime) if startdate else None
    return {'band': validate_band(raw.get('band')), 'stations': stations, 'targets': targets or None,
            'duration': validate_duration(raw.get('duration')),
            'ontarget': _validate_ontarget_fraction(raw.get('ontarget')),
            'startdate': startdate if start_time is not None else None,
            'starttime': starttime if start_time is not None else None, 'start_time': start_time,
            'datarate': validate_datarate(raw.get('datarate')), 'subbands': validate_subbands(raw.get('subbands')),
            'channels': validate_channels(raw.get('channels')),
            'polarizations': validate_polarizations(raw.get('polarizations')),
            'inttime': validate_inttime(raw.get('inttime'))}


def _group_codenames() -> dict[str, set[str]]:
    """Return {group_name: set of station codenames} for grouped stations."""
    groups: dict[str, set[str]] = {}
    for station in observation._STATIONS:
        if station.group:
            groups.setdefault(station.group, set()).add(station.codename)
    return groups


def _passes_validator(value: Any, validator) -> bool:
    """Return True if validator(value) does not raise InvalidObsParams."""
    try:
        validator(value)
        return True
    except InvalidObsParams:
        return False


def validate_url_component(component_id: Any, prop: str, value: Any) -> tuple[bool, Any]:
    """Validate one (component id, property, value) entry from a shared-link ``?config=`` payload.

    Parameters
    ----------
    component_id : str or dict
        Dash component id (string or pattern-matching dict).
    prop : str
        Component property name.
    value : Any
        Raw value from the URL.

    Returns
    -------
    tuple[bool, Any]
        (ok, normalized_value). When ok is False the caller must ignore the value.
    """
    bool_ids = {'switch-band-label', 'switch-specify-epoch', 'switch-specify-e-evn', 'switch-specify-continuum'}
    if isinstance(component_id, dict):
        ctype, index = component_id.get('type'), component_id.get('index')
        if ctype in ('network-switch', 'group-is-selected'):
            return isinstance(value, bool), value
        if ctype == 'group-active-codename':
            return isinstance(value, str) and value in _group_codenames().get(index, set()), value
        return False, None
    if component_id in bool_ids:
        return isinstance(value, bool), value
    if component_id == 'band-slider':
        return _passes_validator(value, band_name_from_index), value
    if component_id == 'duration':
        return _passes_validator(value, validate_duration), value
    if component_id == 'onsourcetime':
        return _passes_validator(value, ontarget_from_percent), value
    if component_id == 'switches-antennas':
        ok = isinstance(value, list) and all(isinstance(c, str) and c in observation._STATIONS for c in value)
        return ok, clean_station_codenames(value) if ok else None
    if component_id == 'store-targets':
        return isinstance(value, list), clean_target_specs(value)
    if component_id == 'startdate':
        if value is None:
            return True, value
        return _passes_validator(value, lambda v: parse_start_time(v, '00:00')), value
    if component_id == 'starttime':
        return _passes_validator(value, lambda v: parse_start_time('2000-01-01', v)), value
    validators = {'datarate': validate_datarate, 'subbands': validate_subbands, 'channels': validate_channels,
                  'pols': validate_polarizations, 'inttime': validate_inttime}
    if component_id in validators:
        return value is not None and _passes_validator(value, validators[component_id]), value
    return False, None
