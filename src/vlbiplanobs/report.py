"""Observation report serialization and file export."""
import json
import logging
import shutil
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Iterable, Optional

import numpy as np
from astropy import units as u
from astropy.time import Time


SUPPORTED_REPORT_EXTENSIONS = {'.pdf', '.txt', '.md', '.json'}

# Outputs shown in the plain-text/Markdown summary (a subset of `observation_report` outputs).
SUMMARY_OUTPUTS = ('frequency', 'wavelength', 'bandwidth', 'on_target_time', 'data_size', 'thermal_noise',
                   'synthesized_beam', 'longest_baseline', 'shortest_baseline', 'bandwidth_smearing',
                   'time_smearing', 'observable_ranges', 'sun_separation')

log = logging.getLogger(__name__)


def validate_report_filename(filename: str) -> Path:
    """Validate a report filename and return it as a path.

    Parameters
    ----------
    filename : str
        Destination filename whose suffix determines the report format.

    Returns
    -------
    Path
        Validated destination path.

    Raises
    ------
    ValueError
        If the filename does not end in .pdf, .txt, .md, or .json.
    """
    path = Path(filename).expanduser()
    if path.suffix.lower() not in SUPPORTED_REPORT_EXTENSIONS:
        allowed = ', '.join(sorted(SUPPORTED_REPORT_EXTENSIONS))
        raise ValueError(f"Unsupported output extension '{path.suffix or '<none>'}'. Use one of: {allowed}")
    return path


def _json_value(value: Any) -> Any:
    """Convert scientific Python values into deterministic JSON-compatible data."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, (float, np.floating)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, u.Quantity):
        raw_value = np.asarray(value.value)
        return {'value': _json_value(raw_value.item() if raw_value.ndim == 0 else raw_value),
                'unit': value.unit.to_string()}
    if isinstance(value, Time):
        isot = value.utc.isot
        return isot if isinstance(isot, str) else np.asarray(isot).tolist()
    if isinstance(value, np.ndarray):
        return [_json_value(item) for item in value.tolist()]
    if isinstance(value, Enum):
        return value.name
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_value(item) for item in value]
    if hasattr(value, 'to_string'):
        return value.to_string()
    return str(value)


def _calculated_value(calculation: Callable[[], Any], errors: dict[str, str], name: str) -> Any:
    """Run one report calculation, recording an explicit error instead of aborting the export."""
    try:
        return _json_value(calculation())
    except Exception as error:
        errors[name] = f"{type(error).__name__}: {error}"
        return None


def observation_report(observation, output_names: Optional[Iterable[str]] = None) -> dict[str, Any]:
    """Return reusable observation inputs and calculated outputs as JSON-compatible data.

    Parameters
    ----------
    observation : VLBIObs
        Computed observation to serialize.
    output_names : iterable of str or None
        Names of the outputs to calculate (e.g. ``SUMMARY_OUTPUTS``). None calculates all of them,
        including the expensive uv values and per-baseline sensitivities.

    Returns
    -------
    dict
        Report with ``inputs``, ``outputs``, and per-calculation ``errors`` sections.

    Raises
    ------
    ValueError
        If ``output_names`` contains an unknown output name.
    """
    scans = {
        block_name: [{
            'name': scan.source.name,
            'coordinates': scan.source.coord.to_string('hmsdms'),
            'type': scan.source.type.name,
            'duration': _json_value(scan.duration),
            'every': scan.every,
        } for scan in block]
        for block_name, block in observation.scans.items()
    }
    inputs = {
        'band': observation.band,
        'stations': list(observation.stations.station_codenames),
        'excluded_stations': dict(observation._excluded_stations),
        'scans': scans,
        'times': _json_value(observation.times),
        'duration': _json_value(observation.duration),
        'datarate': _json_value(observation.datarate),
        'subbands': observation.subbands,
        'channels_per_subband': observation.channels,
        'polarizations': observation.polarizations.value,
        'integration_time': _json_value(observation.inttime),
        'bits_per_sample': _json_value(observation.bitsampling),
        'on_target_fraction': observation.ontarget_fraction,
    }
    errors: dict[str, str] = {}
    calculations = {
        'wavelength': lambda: observation.wavelength,
        'frequency': lambda: observation.frequency,
        'bandwidth': lambda: observation.bandwidth,
        'on_target_time': lambda: observation.ontarget_time,
        'data_size': lambda: observation.datasize(),
        'thermal_noise': lambda: observation.thermal_noise(),
        'synthesized_beam': lambda: observation.synthesized_beam(),
        'longest_baseline': lambda: observation.longest_baseline(),
        'shortest_baseline': lambda: observation.shortest_baseline(),
        'bandwidth_smearing': lambda: observation.bandwidth_smearing(),
        'time_smearing': lambda: observation.time_smearing(),
        'source_observability': lambda: observation.per_source_observable(),
        'block_observability': lambda: observation.can_be_observed(),
        'observable_by_network': lambda: observation.is_observable_by_network(),
        'observable_ranges': lambda: observation.when_is_observable(),
        'sun_separation': lambda: observation.sun_constraint_per_source(),
        'sun_limiting_epochs': lambda: observation.sun_limiting_epochs(),
        'elevations': lambda: observation.per_source_elevations(),
        'uv_values': lambda: observation.get_uv_values(),
        'baseline_sensitivity': lambda: observation.baseline_sensitivity(),
    }
    selected = list(calculations) if output_names is None else list(output_names)
    unknown = [name for name in selected if name not in calculations]
    if unknown:
        raise ValueError(f"Unknown report output(s): {', '.join(unknown)}")
    outputs = {name: _calculated_value(calculations[name], errors, name) for name in selected}
    return {'inputs': _json_value(inputs), 'outputs': outputs, 'errors': errors}


def render_observation_summary(observation, markdown: bool = False) -> str:
    """Render a human-readable plain-text or Markdown observation summary."""
    report = observation_report(observation, output_names=SUMMARY_OUTPUTS)
    inputs = report['inputs']
    outputs = report['outputs']
    heading = '# PlanObs Observation Summary' if markdown else 'PlanObs Observation Summary'
    section = '##' if markdown else ''
    bullet = '- ' if markdown else '  '
    lines = [heading, '', f"{section + ' ' if section else ''}Observation setup"]
    lines.extend([
        f"{bullet}Band: {inputs['band']}",
        f"{bullet}Stations: {', '.join(inputs['stations'])}",
        f"{bullet}Duration: {_display_value(inputs['duration'])}",
        f"{bullet}Times: {_display_value(inputs['times'])}",
        f"{bullet}Data rate: {_display_value(inputs['datarate'])}",
        f"{bullet}Subbands: {inputs['subbands']}",
        f"{bullet}Channels per subband: {inputs['channels_per_subband']}",
        f"{bullet}Polarizations: {inputs['polarizations']}",
        f"{bullet}Integration time: {_display_value(inputs['integration_time'])}",
        f"{bullet}On-target fraction: {inputs['on_target_fraction']}",
        '', f"{section + ' ' if section else ''}Sources",
    ])
    if inputs['scans']:
        for block_name, scans in inputs['scans'].items():
            source_names = ', '.join(f"{scan['name']} ({scan['coordinates']})" for scan in scans)
            lines.append(f"{bullet}{block_name}: {source_names}")
    else:
        lines.append(f"{bullet}No sources defined")
    lines.extend(['', f"{section + ' ' if section else ''}Calculated results"])
    for name in SUMMARY_OUTPUTS:
        lines.append(f"{bullet}{name.replace('_', ' ').title()}: {_display_value(outputs[name])}")
    if report['errors']:
        lines.extend(['', f"{section + ' ' if section else ''}Unavailable calculations"])
        lines.extend(f"{bullet}{name}: {error}" for name, error in report['errors'].items())
    return '\n'.join(lines).rstrip() + '\n'


def _display_value(value: Any) -> str:
    """Format a JSON-compatible report value for human-readable output."""
    if value is None:
        return 'Not available'
    if isinstance(value, dict) and set(value) == {'value', 'unit'}:
        return f"{value['value']} {value['unit']}"
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False)
    return str(value)


def _pdf_observations(observation) -> list:
    """Split a multi-target observation into single-block observations for PDF pages."""
    if len(observation.scans) <= 1:
        return [observation]
    pages = []
    for block_name, block in observation.scans.items():
        page = type(observation)(observation.band, observation.stations, scans={block_name: block},
                                 times=observation.times if observation.fixed_time else None,
                                 duration=observation.duration, datarate=observation.datarate,
                                 subbands=observation.subbands, channels=observation.channels,
                                 polarizations=observation.polarizations.value, inttime=observation.inttime,
                                 ontarget=observation.ontarget_fraction, bits=observation.bitsampling)
        page._excluded_stations = dict(observation._excluded_stations)
        pages.append(page)
    return pages


def write_observation_report(observation, filename: str) -> Path:
    """Write an observation report in the format selected by the filename suffix.

    Parameters
    ----------
    observation : VLBIObs
        Computed observation to export.
    filename : str
        Destination; its suffix (.pdf, .txt, .md, .json) selects the format. Existing files are overwritten.

    Returns
    -------
    Path
        The written file path.

    Raises
    ------
    ValueError
        If the extension is unsupported.
    FileNotFoundError
        If the destination directory does not exist (checked before any slow rendering).
    """
    path = validate_report_filename(filename)
    if not path.parent.is_dir():
        raise FileNotFoundError(f"Output directory '{path.parent}' does not exist.")
    suffix = path.suffix.lower()
    if suffix == '.json':
        content = json.dumps(observation_report(observation), indent=2, ensure_ascii=False) + '\n'
        path.write_text(content, encoding='utf-8')
    elif suffix in {'.txt', '.md'}:
        path.write_text(render_observation_summary(observation, markdown=suffix == '.md'), encoding='utf-8')
    else:
        from vlbiplanobs.gui import outputs
        pages = _pdf_observations(observation)
        try:
            temporary_path = outputs.summary_pdf_for_sources(pages, show_figure=True)
        except Exception as error:
            log.warning("PDF report with figures failed (%s: %s); retrying without figures",
                        type(error).__name__, error)
            temporary_path = outputs.summary_pdf_for_sources(pages, show_figure=False)
        try:
            shutil.move(temporary_path, path)
        except Exception:
            Path(temporary_path).unlink(missing_ok=True)
            raise
    return path
