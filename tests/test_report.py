import json

import pytest
from astropy import units as u

from vlbiplanobs import observation as obs
from vlbiplanobs.cli import VLBIObs
from vlbiplanobs.report import (observation_report, render_observation_summary,
                                validate_report_filename, write_observation_report)


@pytest.fixture
def minimal_observation():
    """Return a deterministic observation that requires no source lookup."""
    stations = obs._STATIONS.filter_antennas(['Ef', 'Wb'])
    return VLBIObs(band='18cm', stations=stations, scans={}, duration=2*u.h,
                   datarate=1024*u.Mbit/u.s, subbands=8, channels=128,
                   polarizations=2, inttime=2*u.s)


def test_validate_report_filename_accepts_only_supported_extensions():
    """Verify report format selection is strict and case-insensitive."""
    for suffix in ('.PDF', '.txt', '.Md', '.JSON'):
        assert validate_report_filename(f'report{suffix}').suffix == suffix
    for filename in ('report', 'report.csv', 'report.pdf.txt.bak'):
        with pytest.raises(ValueError):
            validate_report_filename(filename)


def test_json_report_contains_reusable_inputs_and_calculated_outputs(minimal_observation, tmp_path):
    """Verify JSON export includes setup fields and calculated result fields."""
    output = tmp_path / 'observation.JSON'
    write_observation_report(minimal_observation, str(output))
    report = json.loads(output.read_text(encoding='utf-8'))

    assert report['inputs']['band'] == '18cm'
    assert report['inputs']['stations'] == ['Ef', 'Wb']
    assert report['inputs']['channels_per_subband'] == 128
    assert report['inputs']['subbands'] == 8
    assert report['inputs']['datarate']['value'] == 1024
    assert {'frequency', 'bandwidth', 'thermal_noise', 'data_size', 'elevations'} <= report['outputs'].keys()


def test_text_and_markdown_reports_are_human_readable(minimal_observation, tmp_path):
    """Verify text formats contain the same core observation summary."""
    text_output = tmp_path / 'observation.txt'
    markdown_output = tmp_path / 'observation.md'
    write_observation_report(minimal_observation, str(text_output))
    write_observation_report(minimal_observation, str(markdown_output))

    text = text_output.read_text(encoding='utf-8')
    markdown = markdown_output.read_text(encoding='utf-8')
    for content in (text, markdown):
        assert 'Observation setup' in content
        assert 'Channels per subband: 128' in content
        assert 'Ef, Wb' in content
        assert 'Calculated results' in content
    assert markdown.startswith('# PlanObs Observation Summary')


def test_pdf_report_uses_gui_summary_generator(minimal_observation, tmp_path, monkeypatch):
    """Verify PDF export delegates to the same multi-page generator as the GUI."""
    destination = tmp_path / 'observation.PDF'
    temporary = tmp_path / 'generated.pdf'

    def fake_summary_pdf(observations, show_figure):
        temporary.write_bytes(b'%PDF report')
        return str(temporary)

    from vlbiplanobs.gui import outputs
    monkeypatch.setattr(outputs, 'summary_pdf_for_sources', fake_summary_pdf)

    result = write_observation_report(minimal_observation, str(destination))

    assert result == destination
    assert destination.read_bytes() == b'%PDF report'


def test_rendered_summary_ends_with_newline(minimal_observation):
    """Verify summaries are suitable for normal command-line text files."""
    assert render_observation_summary(minimal_observation).endswith('\n')
    assert observation_report(minimal_observation)['errors'] is not None
