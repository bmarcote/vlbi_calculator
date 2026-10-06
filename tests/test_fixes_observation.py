"""Regression tests for fixes in observation.py and stations.py (shared-state leak, datarate caps,
time-range visibility, empty uv coverage, cache invalidation, key-file injection, catalog reading)."""
import pytest
import numpy as np
from astropy import units as u
from astropy.time import Time
from vlbiplanobs import observation as obs
from vlbiplanobs import stations
from vlbiplanobs import sources


def _obs(codenames: list[str], band: str = '18cm', datarate=1024*u.Mbit/u.s, coords: str = '10h00m00s +40d00m00s',
         fixed: bool = True) -> obs.Observation:
    """Builds an Observation with one target for the given station codenames (module-level shared Stations)."""
    scans = {'target': sources.ScanBlock([sources.Scan(sources.Source('Target', coords))])}
    times = Time('2026-03-01 00:00') + np.arange(0, 24*60, 15)*u.min if fixed else None
    return obs.Observation(band=band, stations=obs._STATIONS.filter_antennas(codenames), scans=scans,
                           times=times, datarate=datarate)


def test_datarate_does_not_mutate_shared_stations():
    """#1: setting the observation datarate must not write onto the shared Station objects."""
    before = {s.codename: s.datarate for s in obs._STATIONS}
    o1 = _obs(['Ef', 'Mc', 'Ys'], datarate=2048*u.Mbit/u.s)
    o2 = _obs(['Ef', 'Mc', 'Ys'], datarate=256*u.Mbit/u.s)
    assert {s.codename: s.datarate for s in obs._STATIONS} == before
    assert o2.station_datarate('Ef') <= 256*u.Mbit/u.s
    assert o1.station_datarate('Ef') > o2.station_datarate('Ef')
    # o1's rms must not change because o2 was created later with a lower datarate
    rms1 = o1.thermal_noise()
    o3 = _obs(['Ef', 'Mc', 'Ys'], datarate=2048*u.Mbit/u.s)
    assert rms1['Target'] == o3.thermal_noise()['Target']


def test_datarate_with_per_network_dict_cap_non_evn():
    """#2: stations with dict max_datarate in a non-EVN network must not raise TypeError."""
    o = obs.Observation(band='6cm', stations=obs._STATIONS.filter_antennas(['Cm', 'Da', 'De']), scans={},
                        datarate=1024*u.Mbit/u.s)
    assert set(o.station_datarates) == {'Cm', 'Da', 'De'}
    assert all(v <= 1024*u.Mbit/u.s for v in o.station_datarates.values())


def test_is_observable_with_times_does_not_touch_cache():
    """#3: is_observable(times=...) with a different length works and doesn't poison the cache."""
    o = _obs(['Ef', 'Mc', 'Ys', 'Hh'])
    other_times = Time('2026-03-01 00:00') + np.arange(0, 6*60, 7)*u.min
    vis_other = o.is_observable(times=other_times)
    assert len(vis_other['target']['Ef']) == len(other_times)
    assert len(o.is_observable()['target']['Ef']) == len(o.times)
    windows = o.when_is_observable(within_time_range=[Time('2026-03-01 00:00'), Time('2026-03-01 12:00')])
    assert 'target' in windows


def test_never_visible_source_beam_and_baselines():
    """#4: a source never visible (dec -85 from EVN) must not crash beam/baseline computations."""
    o = _obs(['Ef', 'Mc', 'Ys', 'Hh', 'Jb2'], coords='10h00m00s -85d00m00s')
    assert 'Target' not in o.synthesized_beam()
    assert 'Target' not in o.longest_baseline()


def test_bitsampling_and_datarate_none_invalidate_rms():
    """#6: changing bitsampling or clearing the datarate invalidates the cached rms."""
    o = _obs(['Ef', 'Mc', 'Ys'])
    rms2 = o.thermal_noise()['Target']
    o.bitsampling = 1*u.bit
    assert o.thermal_noise()['Target'] != rms2
    o.datarate = None
    assert o.thermal_noise() is None
    assert o.station_datarates == {}


def test_wraparound_window_merged_on_reference_day():
    """#5: on the default reference day a window crossing midnight is returned as one window."""
    # RA 0h transits near midnight on the reference day (Sep 21): visible window spans 00:00 UTC
    o = _obs(['Ef', 'Mc', 'Ys', 'Hh'], coords='00h00m00s +20d00m00s', fixed=False)
    windows = o.when_is_observable(min_stations=1)['target']
    assert len(windows) == 1
    assert windows[0][1] - windows[0][0] > 0*u.h
    assert windows[0][1] > o.times[-1]
    o_circ = _obs(['Ef', 'Mc', 'Ys', 'Hh'], coords='12h00m00s +89d00m00s', fixed=False)
    assert len(o_circ.when_is_observable(min_stations=1)['target']) == 1


def test_schedule_file_no_placeholder_injection():
    """#8: user values containing placeholders, quotes or newlines are inserted literally and sanitized."""
    o = _obs(['Ef', 'Mc', 'Ys'])
    key = o.schedule_file(pi_name="O'Brien {SCANS}\nexpcode='HACK'", comments='{SOURCES}')
    assert "piname   = 'OBrien {SCANS} expcode=HACK'" in key
    assert '{SOURCES}' in key
    assert "\nexpcode='HACK'" not in key


def test_custom_network_catalog_is_read(tmp_path):
    """#10a: a custom network catalog file is actually parsed (previously yielded 0 networks)."""
    catalog = tmp_path / 'networks.inp'
    catalog.write_text("[TESTNET]\nname = Test Network\ndefault_antennas = Ef, Mc\n"
                       "max_datarate = 1024\nobserving_bands = 18cm, 6cm\n")
    networks = stations.Stations.get_networks_from_configfile(filename=str(catalog))
    assert list(networks) == ['TESTNET']
    assert networks['TESTNET'].station_codenames == ['Ef', 'Mc']
    with pytest.raises(FileNotFoundError):
        stations.Stations.get_networks_from_configfile(filename=str(tmp_path / 'missing.inp'))


def test_missing_custom_station_catalog_raises(tmp_path):
    """#10b: a missing custom station catalog raises instead of being silently ignored."""
    with pytest.raises(FileNotFoundError):
        stations.Stations(filename=str(tmp_path / 'missing_stations.inp'))


def test_stations_add_and_delitem_deterministic():
    """#10c/d: adding networks keeps a deterministic order; deleting by integer index works."""
    all_st = stations.Stations()
    combined = all_st.filter_antennas(['Ef', 'Mc']) + all_st.filter_antennas(['Mc', 'Ys'])
    assert combined.station_codenames == ['Ef', 'Mc', 'Ys']
    del combined[0]
    assert combined.station_codenames == ['Mc', 'Ys']


def test_is_always_observable_matches_station_constraints():
    """is_always_observable() is derived from is_observable(): it must agree with the astroplan constraints."""
    from vlbiplanobs import cli
    stations = ['Ef', 'Mc', 'Nt', 'Wb', 'Hh', 'Ys', 'T6', 'Ur', 'Sc', 'Mk']
    for target in ('12h29m06.7s +02d03m08.6s', '05h00m00s -60d00m00s', '12h00m00s +88d00m00s'):
        for start_time, duration in ((None, None), (Time('2026-11-02 10:00'), 4 * u.h)):
            obs = cli.main(band='18cm', stations=stations, targets=[target], duration=duration, start_time=start_time,
                           datarate=1024 * u.Mbit / u.s, subbands=8, channels=64, polarizations=4, inttime=2 * u.s,
                           ontarget=0.7)
            always = obs.is_always_observable()
            for blockname, block in obs.scans.items():
                for station in obs.stations:
                    expected = all(bool(station.is_always_observable(obs.times, src)) for src in block.sources())
                    assert always[blockname][station.codename] == expected, (target, station.codename)
                    assert always[blockname][station.codename] == bool(all(obs.is_observable()[blockname]
                                                                           [station.codename]))


def test_sun_separation_is_cached_and_read_only():
    from astropy import coordinates as coord
    times = Time('2026-03-01 00:00') + np.arange(0, 30) * u.day
    src = sources.Source('cache-test', '00h10m00s +01d00m00s')
    expected = src.coord.transform_to(coord.GCRS(obstime=times)).separation(coord.get_sun(times))
    first = src.sun_separation(times)
    assert np.allclose(first.deg, expected.deg, rtol=0, atol=1e-12)
    # same epochs in a new Time object and a new (equal) source: served from the cache
    again = sources.Source('cache-test-2', '00h10m00s +01d00m00s').sun_separation(times.copy())
    assert again is first
    with pytest.raises(ValueError):
        first[0] = 0 * u.deg
    # different epochs or coordinates are computed independently
    assert src.sun_separation(times + 1 * u.day) is not first
    other = sources.Source('cache-test-3', '12h10m00s +01d00m00s').sun_separation(times)
    assert not np.allclose(other.deg, first.deg)
    assert len(src.sun_constraint(15 * u.deg, times=times)) == int(np.sum(expected < 15 * u.deg))
