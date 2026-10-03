"""Regression tests for ObservationScheduler fixes (no obs mutation, real calibrator visibility,
cycle detection, sliding-window placement and SCHED string sanitisation)."""
import numpy as np
import pytest
from astropy import units as u
from astropy.time import Time
from vlbiplanobs import observation as obs
from vlbiplanobs.scheduler import ObservationScheduler, _sched_safe
from vlbiplanobs.sources import Source, SourceType, Scan, ScanBlock


def _make_obs() -> obs.Observation:
    """Build a small phase-referencing observation with eMERLIN stations (triggers the 3C286 block).

    Returns
    -------
    obs.Observation
        Observation with one science block 'B1' over 8 h in 10-min steps.
    """
    ants = ['Ef', 'Jb2', 'Cm', 'O8', 'T6', 'Wb']
    tgt = Source('TGT', '13h20m00s +40d00m00s', source_type=SourceType.TARGET)
    pc = Source('PCAL', '13h22m00s +40d30m00s', source_type=SourceType.PHASECAL)
    scans = {'B1': ScanBlock([Scan(pc, 1.5 * u.min), Scan(tgt, 3.5 * u.min)])}
    return obs.Observation(band='6cm', stations=obs._STATIONS.filter_antennas(ants), scans=scans,
                           times=Time('2024-03-10 18:00', scale='utc') + np.arange(0, 480, 10) * u.min,
                           datarate=1024 * u.Mbit / u.s)


@pytest.fixture(scope='module')
def scheduled():
    """Run the scheduler once on a fresh observation; return (observation, scheduler)."""
    o = _make_obs()
    sch = ObservationScheduler(o, fringefinder_spec=['3C345'])
    sch.schedule()
    return o, sch


def test_schedule_does_not_mutate_observation(scheduled):
    o, sch = scheduled
    assert list(o.scans) == ['B1']
    assert {s.name for s in o.sources()} == {'TGT', 'PCAL'}
    assert any(n.startswith('FF_') for n in sch._scans)
    rms = o.thermal_noise()
    assert rms is not None


def test_calibrator_blocks_use_real_visibility(scheduled):
    o, sch = scheduled
    added = [n for n in sch._scans if n not in o.scans]
    assert added
    for name in added:
        src = sch._scans[name].sources()[0]
        ref = obs.Observation._batch_visibility_erfa(list(o.stations), [src], o.times)[src.name]
        expected = np.sum([ref[st.codename] for st in o.stations], axis=0)
        np.testing.assert_array_equal(sch._vis[name], expected)
        assert not np.allclose(sch._elev[name], 45.0)
    for sb in sch.get_scheduled_blocks():
        if sb.name.startswith(('FF', 'eMERLIN', 'POLCAL')):
            assert sb.n_antennas >= sch.min_ant


def _reference_detect_cycle(scans):
    """Original O(n^2) cycle detection, used as the equivalence oracle."""
    n = len(scans)
    best_clen, best_reps, best_covered = 0, 1, 0
    for clen in range(2, n // 2 + 1):
        pattern = scans[:clen]
        reps, i = 0, 0
        while i + clen <= n:
            if ObservationScheduler._scans_match(pattern, scans[i:i + clen]):
                reps += 1
                i += clen
            else:
                break
        covered = reps * clen
        if reps >= 2 and covered > best_covered:
            best_clen, best_reps, best_covered = clen, reps, covered
    if best_clen > 0:
        return scans[:best_clen], best_reps, scans[best_covered:]
    return scans, 1, []


def test_detect_cycle_matches_reference():
    a, b, c = (Source(n, '10h00m00s +10d00m00s') for n in ('A', 'B', 'C'))
    cases = [
        [Scan(a, 1 * u.min), Scan(b, 3 * u.min)] * 5 + [Scan(a, 1 * u.min)],
        ([Scan(a, 1 * u.min), Scan(b, 3 * u.min), Scan(a, 1 * u.min), Scan(b, 3 * u.min),
          Scan(c, 2 * u.min)] * 3 + [Scan(a, 1 * u.min)]),
        [Scan(a, 60 * u.s), Scan(b, 180.4 * u.s), Scan(a, 60.6 * u.s), Scan(b, 179.5 * u.s), Scan(a, 61.2 * u.s)],
        [Scan(a, 1 * u.min), Scan(b, 1 * u.min), Scan(c, 1 * u.min)],
        [Scan(a, 1 * u.min)],
        [],
    ]
    for scans in cases:
        got = ObservationScheduler._detect_cycle(scans)
        ref = _reference_detect_cycle(scans)
        assert got[1] == ref[1]
        assert [id(s) for s in got[0]] == [id(s) for s in ref[0]]
        assert [id(s) for s in got[2]] == [id(s) for s in ref[2]]


def test_find_best_matches_loop(scheduled):
    _, sch = scheduled
    rng = np.random.default_rng(1)
    name = 'B1'
    vis_orig, elev_orig = sch._vis[name], sch._elev[name]
    try:
        for _ in range(20):
            sch._vis[name] = rng.integers(0, sch._n_ant + 1, sch._n_times)
            elev = rng.uniform(-20, 80, sch._n_times)
            elev[rng.integers(0, sch._n_times, 3)] = np.nan
            sch._elev[name] = elev
            t0, t1 = sch.obs.times[0], sch.obs.times[-1]
            dur = int(rng.integers(5, 90)) * u.min
            got = sch._find_best(name, t0, t1, dur)
            i0, i1 = sch._t2i(t0), sch._t2i(t1)
            steps = max(1, int(np.ceil(dur.to(u.min).value / sch._dt)))
            best_score, best_i = -1.0, -1
            for i in range(i0, i1 - steps + 1):
                n = int(np.min(sch._vis[name][i:i + steps]))
                e = float(np.mean(sch._elev[name][i:i + steps]))
                if n >= sch.min_ant and not np.isnan(e) and n * 100.0 + e > best_score:
                    best_score, best_i = n * 100.0 + e, i
            if best_i < 0:
                assert got is None
            else:
                assert got is not None and sch._t2i(got[0]) == best_i
    finally:
        sch._vis[name], sch._elev[name] = vis_orig, elev_orig


def test_key_file_strips_quotes_and_newlines(scheduled):
    _, sch = scheduled
    assert _sched_safe("O'Brien\nendcover /") == 'OBrien endcover /'
    key = sch.generate_key_file(experiment_code='TEST', pi_name="Evil'\nexpcode='X", comments="a\nendcover /",
                                setup_file='dummy.set')
    assert "piname   = 'Evil expcode=X'" in key
    assert '\nendcover /\n' in key and key.count('endcover /') == 2


def _scheduler_for(ants: list[str]) -> ObservationScheduler:
    """Build an ObservationScheduler (not run) for the given antennas on a simple single-target observation."""
    tgt = Source('TGT', '13h20m00s +40d00m00s', source_type=SourceType.TARGET)
    o = obs.Observation(band='18cm', stations=obs._STATIONS.filter_antennas(ants),
                        scans={'B1': ScanBlock([Scan(tgt, 5 * u.min)])},
                        times=Time('2024-03-10 18:00', scale='utc') + np.arange(0, 480, 10) * u.min,
                        datarate=1024 * u.Mbit / u.s)
    return ObservationScheduler(o)


def test_emerlin_3c286_not_triggered_by_jb2_alone():
    """Jodrell Bank (Jb1/Jb2) observes regularly in the EVN: on its own it does not trigger the 3C286 scan."""
    assert not _scheduler_for(['Ef', 'Jb2', 'O8', 'T6', 'Wb', 'Mc'])._has_emerlin()
    assert not _scheduler_for(['Ef', 'Jb1', 'O8', 'T6', 'Wb', 'Mc'])._has_emerlin()
    assert _scheduler_for(['Ef', 'Jb2', 'Cm', 'O8', 'T6', 'Wb'])._has_emerlin()


def test_evn_only_schedule_has_no_3c286():
    """End to end: scheduling an EVN-only array (with Jb2) never adds the eMERLIN_3C286 block."""
    sch = _scheduler_for(['Ef', 'Jb2', 'O8', 'T6', 'Wb', 'Mc'])
    sch.schedule()
    assert 'eMERLIN_3C286' not in sch._scans
