import pytest
import numpy as np
# from astropy import coordinates as coord
from astropy import units as u
from astropy.time import Time
# from rich import print as rprint
from vlbiplanobs import observation as obs
from vlbiplanobs import stations
# from vlbiplanobs import sources
from vlbiplanobs import cli


def test_polarizations():
    assert obs.Polarization.SINGLE == obs.Polarization(1)
    assert obs.Polarization.SINGLE == obs.Polarization['SINGLE']
    assert obs.Polarization.DUAL == obs.Polarization(2)
    assert obs.Polarization.DUAL == obs.Polarization['DUAL']
    assert obs.Polarization.FULL == obs.Polarization(4)
    assert obs.Polarization.FULL == obs.Polarization['FULL']


def test_observation_init():
    with pytest.raises(TypeError):
        o = obs.Observation()

    o = obs.Observation(band='18cm', stations=stations.Stations(), scans={})
    o.polarizations = 'dual'
    assert o.polarizations.value == 2
    o.polarizations = 4
    assert o.polarizations == obs.Polarization.FULL
    evn6 = ['Ef', 'Jb2', 'O8', 'Hh', 'T6', 'Wb', 'Sv', 'Zc', 'Pa', 'Mp', 'Ho', 'Nl', 'Pt', 'Sc', 'Kp', 'Hn']

    for tar_name, tar_coord in (('Target', '10h58m29.6s +81d33m58.8s'),
                                (['Tar1', 'Tar2'],
                                 ['10h58m29.6s +81d33m58.8s', '20h58m29.6s +71d33m58.8s'])):
        if isinstance(tar_name, list):
            scans = [obs.Scan(source=obs.Source(aname, acoord))
                     for aname, acoord in zip(tar_name, tar_coord)]
            target = [s.source for s in scans]
        else:
            target = obs.Source(tar_name, tar_coord)
            scans = [obs.Scan(source=target)]

        o = obs.Observation(band='18cm', stations=obs._STATIONS.filter_antennas(evn6),
                            scans={'block1': obs.ScanBlock(scans)})
        o.times = Time('2020-06-15 20:00', scale='utc') + np.arange(0, 720, 10)*u.min
        assert o.band == '18cm'
        o.datarate = 1024*u.Mbit/u.s
        o.subbands = 8
        o.channels = 32
        o.polarizations = 2
        o.polarizations = 'full'
        with pytest.raises(ValueError):
            o.polarizations = 0
            o.polarizations = -5
            o.polarizations = 'efww'
            o.polarizations = ''
            o.polarizations = None

        o.inttime = 2*u.s

        if isinstance(tar_name, list):
            assert len(o.sources()) == len(target)
        else:
            assert len(o.sources()) == 1

        # assert len(o.sources) == len(tar_name) if isinstance(tar_name, list) else 1
        # all_stations = stations.Network.get_stations_from_configfile()
        # o.stations_from_codenames(evn6)
        _ = o.elevations()
        temp = o.is_observable()
        assert len(temp) == 1 and len(list(temp.values())[0]) == len(evn6)
        # temp = obs.Observation.guest_times_for_source(target if not isinstance(target, list) else target[0],
        #                                               obs.Network.get_stations_from_configfile(codenames=evn6))
        # beam = o.synthesized_beam)
        temp = o.is_always_observable()
        assert len(temp) == 1
        assert list(temp.values())[0]['Ef']
        # _ = o.thermal_noise()
        # _ = o.get_uv_array()
        # _ = o.get_dirtymap(robust='natural')
        # _ = o.get_dirtymap(robust='uniform')
        # fig, ax = plt.subplots()

        assert 'Ef' in o.stations.station_codenames, f"Stations found: {o.stations.station_codenames}"
        assert 'Hh' in o.stations.station_codenames, f"Stations found: {o.stations.station_codenames}"
        assert 'Jb2' in o.stations.station_codenames, f"Stations found: {o.stations.station_codenames}"
        assert 'Nl' in o.stations.station_codenames, f"Stations found: {o.stations.station_codenames}"
        assert 'Pa' in o.stations.station_codenames, f"Stations found: {o.stations.station_codenames}"
        with pytest.raises(AssertionError):
            assert 'Cm' in o.stations.station_codenames, f"Stations found: {o.stations.station_codenames}"
            assert 'Ku' in o.stations.station_codenames, f"Stations found: {o.stations.station_codenames}"


# TODO:
# test uv coverage: same number of points in +x,+y than -x,-y., and the other quarter.
# test rms: scaling by time as sqrt(), same for banddith, etc.

def test_thermal_noise():
    evn6 = ['Ef', 'Jb2', 'O8', 'T6', 'Wb']
    o = cli.VLBIObs(band='18cm', stations=obs._STATIONS.filter_antennas(evn6), duration=10*u.h, scans={},
                    datarate=1024*u.Mbit/u.s)
    rmss = []
    durations = np.array([1, 5, 10, 20])*u.h
    for durs in durations:
        o.duration = durs
        rmss.append(o.thermal_noise())

    assert len(rmss) == len(durations), f"Something went wrong here {rmss=}"

    rmss = np.array([rms.value for rms in rmss])
    assert all([(rmss[i+1] - rms) < 0.0 for i, rms in enumerate(rmss[:-1])]), \
        "Larger obs time should produce lower rms."
    assert all(np.abs(rmss[1:] - np.array([rms/np.sqrt(t-durations[0]).value
                      for t, rms in zip(durations[1:], rmss[1:])])) < 1e-5), \
        f"The returned rms are (for {durations}): {rmss} Jy/beam."

    evn6 = ['Ef', 'Jb2', 'T6', 'Wb']
    o = cli.VLBIObs(band='18cm', stations=obs._STATIONS.filter_antennas(evn6), duration=10*u.h, scans={},
                    datarate=1024*u.Mbit/u.s)
    assert o.thermal_noise() > rmss[2]*u.Jy/u.beam
