"""Tests for the Network Monitoring Experiment (NME) planning in vlbiplanobs.nme."""
import pytest
from astropy import units as u
from astropy.time import Time
from vlbiplanobs import nme
from vlbiplanobs import observation as obs


def _check_contiguous(scans, duration_s):
    assert scans[0].slot_start_s == 0 and scans[0].gap_s == 0
    assert scans[-1].stop_s == duration_s
    assert all(a.stop_s == b.slot_start_s for a, b in zip(scans, scans[1:]))
    assert all(s.dur_s > 0 for s in scans)


def test_grab_offsets_long_observation_every_30_min():
    assert nme.grab_offsets(3 * 3600) == [600, 2400, 4200, 6000, 7800, 9600]


def test_grab_offsets_short_observation_every_15_min_plus_final():
    assert nme.grab_offsets(3600) == [300, 1200, 2100, 3000, 3600 - nme.GRAB_TAIL_S]


def test_grab_offsets_too_short():
    with pytest.raises(ValueError):
        nme.grab_offsets(10 * 60)


@pytest.mark.parametrize('hours', [0.5, 1, 2, 2.5, 3, 7.3, 12])
def test_plan_scans_covers_full_time_and_aligns_grabs(hours):
    duration_s = int(hours * 3600)
    grabs = nme.grab_offsets(duration_s)
    scans = nme.plan_scans(duration_s, grabs)
    _check_contiguous(scans, duration_s)
    grab_scans = [s for s in scans if s.grab_s is not None]
    assert [s.grab_s for s in grab_scans] == grabs
    for s in grab_scans:
        assert s.stop_s == s.grab_s + nme.GRAB_TAIL_S
        assert s.dur_s >= nme.GRAB_TAIL_S + nme.GRAB_DATA_S


def test_plan_scans_matches_reference_3h_layout():
    scans = nme.plan_scans(3 * 3600, nme.grab_offsets(3 * 3600))
    # After the first grab, scans follow the N25L1 pattern: gap=5:00 dur=10:00, last one dur=13:00.
    assert [(s.gap_s, s.dur_s) for s in scans[2:-1]] == [(300, 600)] * 10
    assert (scans[-1].gap_s, scans[-1].dur_s) == (300, 780)


def test_plan_nme_and_key_file():
    stations = obs._STATIONS.filter_antennas(['Ef', 'Wb', 'O8', 'Tr', 'Nt'])
    start = Time('2025-02-20 12:00', scale='utc')
    scans, sources, counts, times = nme.plan_nme(stations, start, 3 * u.hour, '18cm')
    assert all(s.source is not None and s.n_visible == len(stations) for s in scans)
    key = nme.generate_nme_key_file(scans, stations, start, '18cm', 'n25l1', "'evn18cm.set'",
                                    datarate_mbps=1024)
    assert "expcode  = 'N25L1'" in key
    assert key.count("grabto='FILE' grabtime=2,118") == 6
    assert "grabto='NONE'" in key
    assert 'start = 12:00:00' in key and 'month = 02' in key
    assert '{' not in key.split('srccat /')[1]
    assert '12:10:00 (scan  2, 2 sec,' in key


def test_resolve_fringe_finders_from_rfc():
    sources = nme.resolve_fringe_finders(['J0237+2848', '3C84'])
    assert [s.name for s in sources] == ['J0237+2848', 'J0319+4130']
