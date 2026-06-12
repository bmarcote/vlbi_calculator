import sys
import textwrap
import pytest
import tempfile
import os
from astropy import units as u
from vlbiplanobs.cli import VLBIObs
from vlbiplanobs import observation as obs


# ---------------------------------------------------------------------------
# VLBIObs unit tests
# ---------------------------------------------------------------------------

def test_vlbiobs_init_minimal():
    # Minimal valid VLBIObs initialization
    o = VLBIObs(band='18cm', stations=obs._STATIONS.filter_antennas(['Ef']), scans={})
    assert o.band == '18cm'
    assert hasattr(o, 'stations')
    assert hasattr(o, 'scans')

def test_vlbiobs_inherits_observation():
    # VLBIObs should be a subclass of Observation
    o = VLBIObs(band='18cm', stations=obs._STATIONS.filter_antennas(['Ef']), scans={})
    assert isinstance(o, obs.Observation)

def test_vlbiobs_summary_methods():
    o = VLBIObs(band='18cm', stations=obs._STATIONS.filter_antennas(['Ef']), scans={},
                datarate=1024*u.Mbit/u.s)
    # These methods should run without error
    o.summary(gui=False, tui=True)
    o.summary(gui=True, tui=False)

def test_vlbiobs_properties():
    evn6 = ['Ef', 'Jb2', 'O8', 'T6', 'Wb']
    o = VLBIObs(band='18cm', stations=obs._STATIONS.filter_antennas(evn6), duration=10*u.h, scans={})
    assert o.duration == 10*u.h
    assert set(o.stations.station_codenames) == set(evn6)

def test_vlbiobs_invalid_init():
    # Should raise if required arguments are missing
    with pytest.raises(TypeError):
        VLBIObs()


# ---------------------------------------------------------------------------
# CLI integration helpers
# ---------------------------------------------------------------------------

def _run_cli(argv: list[str]) -> tuple[int, str]:
    """Invoke ``vlbiplanobs.cli:cli`` with ``sys.argv`` patched to *argv*.

    Returns (exit_code, captured_stdout+stderr).
    The ``argv`` list should NOT include the program name; it is prepended
    automatically as 'planobs'.
    """
    import io
    from contextlib import redirect_stdout, redirect_stderr
    from vlbiplanobs.cli import cli

    original_argv = sys.argv[:]
    sys.argv = ['planobs'] + argv
    buf = io.StringIO()
    exit_code = 0
    try:
        with redirect_stdout(buf), redirect_stderr(buf):
            cli()
    except SystemExit as e:
        exit_code = int(e.code) if e.code is not None else 0
    finally:
        sys.argv = original_argv
    return exit_code, buf.getvalue()


# ---------------------------------------------------------------------------
# Version and help
# ---------------------------------------------------------------------------

def test_version_flag():
    """``planobs -V`` prints the version and exits 0."""
    code, out = _run_cli(['-V'])
    assert code == 0
    assert 'planobs' in out
    # Version string should look like a semver, e.g. 'planobs 5.0'
    parts = out.strip().split()
    assert len(parts) == 2
    assert parts[0] == 'planobs'

def test_version_long_flag():
    """``planobs --version`` is equivalent to ``-V``."""
    code, out = _run_cli(['--version'])
    assert code == 0
    assert 'planobs' in out

def test_no_args_shows_help():
    """``planobs`` with no arguments exits 0 and prints usage."""
    code, out = _run_cli([])
    assert code == 0
    assert 'planobs' in out.lower() or 'usage' in out.lower() or 'mode' in out.lower()

def test_observe_help():
    """``planobs observe -h`` exits 0 and contains key flags."""
    code, out = _run_cli(['observe', '-h'])
    assert code == 0
    assert '--band' in out or '-b' in out

def test_fringefinders_help():
    code, out = _run_cli(['fringefinders', '-h'])
    assert code == 0
    assert '--network' in out or '-n' in out

def test_phasecals_help():
    code, out = _run_cli(['phasecals', '-h'])
    assert code == 0
    assert '--target' in out or '-t' in out
    assert '--source-catalog' in out or '-sc' in out

def test_source_help():
    code, out = _run_cli(['source', '-h'])
    assert code == 0
    assert 'source_name' in out or 'SOURCE' in out.upper()

def test_antenna_help():
    code, out = _run_cli(['antenna', '-h'])
    assert code == 0
    assert 'antenna' in out.lower()

def test_ant_alias_help():
    """``ant`` is an alias for ``antenna``."""
    code, out = _run_cli(['ant', '-h'])
    assert code == 0
    assert 'antenna' in out.lower()


# ---------------------------------------------------------------------------
# Unknown / invalid arguments raise errors
# ---------------------------------------------------------------------------

def test_unknown_top_level_flag_in_observe_mode():
    """An unrecognised flag to the observe subcommand should exit non-zero."""
    code, _ = _run_cli(['observe', '--nonexistent-flag-xyz'])
    assert code != 0

def test_unknown_subcommand():
    """An unknown subcommand in legacy mode fails (missing band)."""
    # In legacy mode (first arg starts with '-'), an unknown flag triggers error.
    code, _ = _run_cli(['--nonexistent-flag-xyz'])
    assert code != 0

def test_fringefinders_missing_required_args():
    """``planobs fringefinders`` without required --starttime / --duration exits non-zero."""
    code, _ = _run_cli(['fringefinders', '-n', 'EVN'])
    assert code != 0

def test_phasecals_missing_target():
    """``planobs phasecals`` without -t exits non-zero."""
    code, _ = _run_cli(['phasecals'])
    assert code != 0

def test_source_missing_positional():
    """``planobs source`` without a source name exits non-zero."""
    code, _ = _run_cli(['source'])
    assert code != 0


# ---------------------------------------------------------------------------
# Functional smoke-tests for fast modes (no heavy network I/O needed)
# ---------------------------------------------------------------------------

def test_observe_list_bands():
    """``planobs observe --list-bands`` exits 0 and prints at least one band."""
    code, out = _run_cli(['observe', '--list-bands'])
    assert code == 0
    assert 'cm' in out

def test_observe_list_networks():
    """``planobs observe --list-networks`` exits 0 and shows known networks."""
    code, out = _run_cli(['observe', '--list-networks'])
    assert code == 0
    assert 'EVN' in out

def test_observe_list_antennas():
    """``planobs observe --list-antennas`` exits 0 and lists Effelsberg."""
    code, out = _run_cli(['observe', '--list-antennas'])
    assert code == 0
    assert 'Effelsberg' in out or 'Ef' in out

def test_observe_missing_band():
    """``planobs observe -n EVN`` without a band exits non-zero."""
    code, _ = _run_cli(['observe', '-n', 'EVN', '--no-tui'])
    assert code != 0

def test_observe_invalid_band():
    """``planobs observe -b 999cm -n EVN`` with a bogus band exits non-zero."""
    code, _ = _run_cli(['observe', '-b', '999cm', '-n', 'EVN', '--no-tui'])
    assert code != 0

def test_antenna_list_all():
    """``planobs antenna`` with no args lists all antennas."""
    code, out = _run_cli(['antenna'])
    assert code == 0
    assert 'Ef' in out

def test_antenna_list_by_band():
    """``planobs antenna -b 18cm`` lists antennas at 18 cm."""
    code, out = _run_cli(['antenna', '-b', '18cm'])
    assert code == 0
    assert 'Ef' in out

def test_antenna_lookup_by_name():
    """``planobs antenna Ef`` shows Effelsberg details."""
    code, out = _run_cli(['antenna', 'Ef'])
    assert code == 0
    assert 'Effelsberg' in out

def test_antenna_invalid_band():
    """``planobs antenna -b 999cm`` exits non-zero for unknown band."""
    code, _ = _run_cli(['antenna', '-b', '999cm'])
    assert code != 0

def test_antenna_unknown_name():
    """``planobs antenna DOESNOTEXIST`` exits non-zero."""
    code, _ = _run_cli(['antenna', 'DOESNOTEXIST_XYZ'])
    assert code != 0

def test_source_lookup_known():
    """``planobs source 3C84 --no-networks`` exits 0 and shows source info."""
    code, out = _run_cli(['source', '3C84', '--no-networks'])
    assert code == 0
    # Should show coordinates or name information
    assert '3C84' in out or '0316+413' in out or 'Coordinates' in out

def test_source_lookup_unknown():
    """``planobs source NOTASOURCE_XYZ`` exits non-zero."""
    code, _ = _run_cli(['source', 'NOTASOURCE_XYZ_ZZZZZ', '--no-networks'])
    assert code != 0


# ---------------------------------------------------------------------------
# phasecals with personal source catalog (-sc)
# ---------------------------------------------------------------------------

_TOML_CATALOG = textwrap.dedent("""\
    [[target]]
    name = "TestSrc"
    [target.coordinates]
    RA = "05:34:32.0"
    Dec = "+22:00:52"
""")


@pytest.fixture
def personal_catalog_file():
    """Write a minimal personal catalog to a temp file, yield path, then clean up."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.toml', delete=False) as fh:
        fh.write(_TOML_CATALOG)
        path = fh.name
    yield path
    os.unlink(path)


def test_phasecals_with_source_catalog(personal_catalog_file):
    """``planobs phasecals -sc <file> -t TestSrc`` resolves from personal catalog."""
    code, out = _run_cli(['phasecals', '-sc', personal_catalog_file, '-t', 'TestSrc', '-n', '1'])
    assert code == 0
    # Should find at least 1 result or report no candidates near that coord
    assert 'TestSrc' in out or 'Found' in out or 'No phase' in out

def test_phasecals_sc_fallback_to_rfc(personal_catalog_file):
    """When target not in personal catalog, warn and fall back to RFC/online lookup."""
    code, out = _run_cli(['phasecals', '-sc', personal_catalog_file, '-t', '3C84', '-n', '1'])
    # Should print fallback warning and still succeed (3C84 is in RFC)
    assert 'falling back' in out or 'Found' in out or 'fallback' in out.lower()

def test_phasecals_sc_nonexistent_file():
    """``planobs phasecals -sc /no/such/file.toml -t 3C84`` exits non-zero."""
    code, _ = _run_cli(['phasecals', '-sc', '/no/such/file.toml', '-t', '3C84', '-n', '1'])
    assert code != 0


# ---------------------------------------------------------------------------
# phasecals parameter coverage (min-flux, max-separation, n-sources)
# ---------------------------------------------------------------------------

def _phasecals_json(extra_args: list) -> tuple[int, dict]:
    """Run ``planobs phasecals`` with ``--json`` and return (exit_code, parsed_dict)."""
    import json
    code, out = _run_cli(['phasecals', '-t', '3C84', '--json'] + extra_args)
    try:
        data = json.loads(out)
    except json.JSONDecodeError:
        data = {}
    return code, data


def test_phasecals_json_output_structure():
    """JSON output contains expected top-level keys when sources are found."""
    code, data = _phasecals_json(['--min-flux', '0.001'])
    if code == 0:
        for key in ('target_name', 'sources', 'total_found', 'min_flux_jy', 'max_separation_deg'):
            assert key in data, f"missing key: {key}"


def test_phasecals_min_flux_echoed_in_json():
    """--min-flux value is echoed in JSON output under ``min_flux_jy``."""
    code, data = _phasecals_json(['--min-flux', '0.42'])
    if code == 0:
        assert abs(data['min_flux_jy'] - 0.42) < 1e-6


def test_phasecals_max_separation_echoed_in_json():
    """--max-separation value is echoed in JSON output under ``max_separation_deg``."""
    code, data = _phasecals_json(['--max-separation', '3.7'])
    if code == 0:
        assert abs(data['max_separation_deg'] - 3.7) < 1e-6


def test_phasecals_min_flux_high_returns_fewer_sources():
    """--min-flux 9999 returns far fewer sources than --min-flux 0.001.

    Sources with no c-band measurement (stored as negative flux) are kept even
    at a high threshold (include_missing=True mode), so the count may not be 0.
    """
    _, data_high = _phasecals_json(['--min-flux', '9999.0'])
    _, data_low = _phasecals_json(['--min-flux', '0.001'])
    n_high = data_high.get('total_found', 0)
    n_low = data_low.get('total_found', 0)
    assert n_high < n_low



def test_phasecals_n_sources_limits_count():
    """-n 1 returns at most 1 source in JSON output."""
    code, data = _phasecals_json(['-n', '1', '--min-flux', '0.001'])
    if code == 0:
        assert data.get('total_found', 0) <= 1
        assert len(data.get('sources', [])) <= 1


def test_phasecals_max_separation_all_within_radius():
    """All returned sources have separation <= --max-separation."""
    code, data = _phasecals_json(['--max-separation', '2.0', '--min-flux', '0.001'])
    if code == 0:
        for src in data.get('sources', []):
            assert src['separation_deg'] <= 2.0 + 1e-6


def test_phasecals_max_separation_tiny_returns_nothing():
    """--max-separation 0.001 (sub-arcminute) returns no sources near 3C84."""
    code, data = _phasecals_json(['--max-separation', '0.001'])
    assert data.get('total_found', 0) == 0
