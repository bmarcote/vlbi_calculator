"""Regression tests for CLI fixes: time grid, source-resolution errors, rich markup escaping,
None data rates, and SCHED experiment code / string validation."""
import argparse
import copy

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates.name_resolve import NameResolveError
from astropy.time import Time

from vlbiplanobs import cli, nme, sources
from vlbiplanobs import observation as obs
from vlbiplanobs.cli import VLBIObs

MARKUP_BREAKER = '[/x]'
START = Time('2025-03-15 08:00', scale='utc')


def _fake_source(name: str = 'FAKE') -> sources.Source:
    """Return an offline Source (no name resolution)."""
    return sources.Source(name=name, coordinates=SkyCoord('12h00m00s +30d00m00s'))


class _EmptyCatalog:
    """Stand-in RFC catalog where no source is found."""

    def get_source(self, name):
        """Return None for every name."""
        return None


# 1. Time grid and duration bounds -------------------------------------------------------------
@pytest.mark.parametrize('duration', [2*u.h, 96*u.min, 7*u.min])
def test_main_time_grid_is_strictly_increasing_and_ends_at_duration(duration):
    """The sampling grid has no duplicates, is monotonic and ends exactly at start + duration."""
    o = cli.main(band='18cm', stations=['Ef', 'Wb'], start_time=START, duration=duration,
                 datarate=1024*u.Mbit/u.s)
    offsets = (o.times - START).to(u.min).value
    assert np.all(np.diff(offsets) > 0)
    assert offsets[0] == pytest.approx(0.0)
    assert offsets[-1] == pytest.approx(duration.to(u.min).value)


@pytest.mark.parametrize('duration', [0*u.h, -1*u.h, 97*u.h])
def test_main_rejects_out_of_bounds_duration(duration):
    """Durations <= 0 or > 96 h raise ValueError before any computation."""
    with pytest.raises(ValueError, match='duration'):
        cli.main(band='18cm', stations=['Ef', 'Wb'], duration=duration, datarate=1024*u.Mbit/u.s)


# 2. Source-resolution errors ------------------------------------------------------------------
def test_resolve_calibrators_skips_name_resolve_error(monkeypatch):
    """A NameResolveError from online lookup is reported and the name skipped, not raised."""
    def raise_nre(*args, **kwargs):
        raise NameResolveError('offline')
    monkeypatch.setattr(sources.Source, 'source_from_str', staticmethod(raise_nre))
    resolved = cli._resolve_calibrators(['NOPE', '[x]NOPE'], _fake_source(), '6cm',
                                        sources.SourceType.PHASECAL, catalog=_EmptyCatalog())
    assert resolved == []


def test_resolve_calibrators_markup_name_gives_value_error_not_markup_error():
    """'[/x]' is parsed as 'name/coordinates'; the bad coordinates raise a plain ValueError."""
    with pytest.raises(ValueError, match='coordinate'):
        cli._resolve_calibrators([MARKUP_BREAKER], _fake_source(), '6cm', sources.SourceType.PHASECAL,
                                 catalog=_EmptyCatalog())


def test_main_wraps_only_target_resolution_errors(monkeypatch):
    """Target resolution failures become ValueError with the cause chained."""
    def raise_nre(*args, **kwargs):
        raise NameResolveError('offline')
    monkeypatch.setattr(sources.Source, 'source_from_str', staticmethod(raise_nre))
    with pytest.raises(ValueError, match="Source 'NOPE' not found") as excinfo:
        cli.main(band='6cm', stations=['Ef', 'Wb'], targets=['NOPE'], datarate=1024*u.Mbit/u.s)
    assert isinstance(excinfo.value.__cause__, NameResolveError)


def test_main_does_not_rewrite_unrelated_errors(monkeypatch):
    """Errors raised after the target is resolved are not disguised as 'source not found'."""
    monkeypatch.setattr(sources.Source, 'source_from_str', staticmethod(lambda *a, **k: _fake_source()))

    def boom(*args, **kwargs):
        raise RuntimeError('calibrator failure')
    monkeypatch.setattr(cli, '_resolve_calibrators', boom)
    monkeypatch.setattr(cli, '_rfc_catalog_for_band', lambda band: _EmptyCatalog())
    with pytest.raises(RuntimeError, match='calibrator failure'):
        cli.main(band='6cm', stations=['Ef', 'Wb'], targets=['FAKE'], phasecal_names=['X'],
                 datarate=1024*u.Mbit/u.s)


# 3. Rich markup escaping ----------------------------------------------------------------------
def test_get_stations_unknown_station_with_markup_raises_value_error():
    """A station name that looks like a closing markup tag must not crash rich."""
    with pytest.raises(ValueError, match='not known'):
        cli.get_stations('18cm', list_stations=[MARKUP_BREAKER])


def test_antenna_command_with_markup_name_exits_cleanly():
    """'planobs antenna [/x]' reports 'not found' and exits with status 1."""
    args = argparse.Namespace(antenna_name=MARKUP_BREAKER, band=None)
    with pytest.raises(SystemExit) as excinfo:
        cli.handle_antenna_command(args)
    assert excinfo.value.code == 1


def test_summary_with_markup_block_and_source_names():
    """Block and source names with markup-like text are printed literally."""
    src = _fake_source(MARKUP_BREAKER)
    o = VLBIObs(band='18cm', stations=obs._STATIONS.filter_antennas(['Ef', 'Wb']),
                scans={MARKUP_BREAKER: sources.ScanBlock([sources.Scan(src, duration=5*u.min)])},
                datarate=1024*u.Mbit/u.s)
    o.summary(gui=False, tui=True)


# 4. None data rates ---------------------------------------------------------------------------
def test_summary_without_observation_datarate():
    """An observation without data rate prints its summary without TypeError."""
    o = VLBIObs(band='18cm', stations=obs._STATIONS.filter_antennas(['Ef', 'Wb']), scans={})
    o.summary(gui=False, tui=True)


def test_summary_with_station_without_datarate():
    """A station with an undefined data rate is not compared against the observation data rate."""
    stations = copy.deepcopy(obs._STATIONS.filter_antennas(['Ef', 'Wb']))
    o = VLBIObs(band='18cm', stations=stations, scans={}, datarate=1024*u.Mbit/u.s)
    o.stations['Ef'].datarate = None
    o.summary(gui=False, tui=True)


# 6. SCHED experiment code and quoted strings --------------------------------------------------
@pytest.mark.parametrize('sched, expected', [('em179a', ('em179a.key', 'EM179A')),
                                             ('out/n25l1.key', ('out/n25l1.key', 'N25L1')),
                                             ('a.b', ('a.b.key', 'A.B'))])
def test_sched_paths(sched, expected):
    """Valid --sched values give the .key filename and uppercased code; invalid ones raise."""
    if '.' in expected[1]:
        with pytest.raises(ValueError, match='Invalid experiment code'):
            cli._sched_paths(sched)
    else:
        assert cli._sched_paths(sched) == expected


@pytest.mark.parametrize('sched', ["x'y", 'bad-name', 'with space', '.key'])
def test_sched_paths_rejects_invalid_codes(sched):
    """Experiment codes with characters outside [A-Za-z0-9_] are rejected."""
    with pytest.raises(ValueError, match='Invalid experiment code'):
        cli._sched_paths(sched)


@pytest.mark.parametrize('bad_name', ["O'Brien", 'two\nlines'])
def test_nme_rejects_unsafe_source_names(bad_name):
    """Source names with quotes or newlines cannot be written into SCHED quoted fields."""
    src = _fake_source(bad_name)
    with pytest.raises(ValueError, match='SCHED'):
        nme._source_catalog_line(src)
    scan = nme.NMEScan(slot_start_s=0, gap_s=0, stop_s=900, source=src, n_visible=2)
    with pytest.raises(ValueError, match='SCHED'):
        nme.scan_lines([scan], START, 2)


def test_nme_accepts_regular_source_name():
    """A normal source name is written unchanged."""
    line = nme._source_catalog_line(_fake_source('J1200+3000'))
    assert line.startswith("source='J1200+3000'")
