import pytest
import numpy as np
from vlbiplanobs import observation as obs
# from vlbiplanobs import stations
from astropy import coordinates as coord
from astropy import units as u
from astropy.time import Time
from rich import print as rprint


def test_source():
    """Tests the Source class
    """
    with pytest.raises(AttributeError):
        obs.Source()

    with pytest.raises(ValueError):
        obs.Source('a_name', 'coordinates')
        obs.Source('a_name', '00h00m00s xxhqqmees')
        obs.Source('a_name', '40h00m00s 10d10m10s')
        obs.Source('a_name', '00:00:00 100:00:00')
        obs.Source('a_name', '-30:00:00 10:00:00')
        # obs.Source('Cyg X-1')
        obs.Source(name='WQQWE')
        obs.Source('WQQWE')

    s1 = obs.Source('a_name', '10h20m10s 40d30m10s')
    s2 = obs.Source('a_name', '10:20:10 40:30:10', unit=(u.hourangle, u.deg))
    try:
        s3 = obs.Source(name='Cyg X-1')
        assert s3.name == 'Cyg X-1'
        s3.coord.to_string('hmsdms')
    except coord.name_resolve.NameResolveError:
        rprint("[yellow]WARNING: seems like no internet connection is available. " \
               "Unable to resolve source names[/yellow]")

    assert s1.name == 'a_name'
    assert s2.name == 'a_name'
    assert s1.coord == s2.coord

    # checking that the distance to the Sun is right
    # times are spring equinox, summer solstice fall equinox
    times = Time(['2024-03-19 22:06', '2024-06-20 15:51', '2024-09-22 07:44'], scale='utc')
    s1 = obs.Source('equinox source', '0h0m0s 0d0m0s')
    # print(f"\n\nSeparations: {s1.sun_separation(times)}\n\n")
    assert all(np.abs(s1.sun_separation(times) - np.array([0.0, 90.0, 180])*u.deg) < 1*u.deg)
    # print(f"\n\nConstraints: {s1.sun_constraint(20*u.deg, times=times)}\n\n")
    assert len(s1.sun_constraint(20*u.deg, times=times)) == 1
    s1 = obs.Source('polar_source', '0h0m0s 90d0m0s')
    # print(f"\n\nSeparations: {s1.sun_separation(times)}\n\n")
    assert all(np.abs(s1.sun_separation(times[(0, 2),]) - np.array([90.0, 90.0])*u.deg) < 1*u.deg)

    # Playing with multiple sources
    s = obs.Source(['s2', 's3'], ['10h20m10s 40d30m10s', '20h30m0s 20d20m10s'])
    sep = s.sun_separation(times)
    assert len(sep) == 2
    assert (len(sep[0]) == len(times)) and (len(sep[1]) == len(times))
    s.sun_constraint(20*u.deg)




def test_source_from_names():
    try:
        assert isinstance(obs.Source.get_source_from_name('Cyg X-1'), obs.Source)
        assert isinstance(obs.Source.get_coordinates_from_name('Cyg X-1'), coord.SkyCoord)

        src_name1 = '2358+390'
        src_name2 = 'J0000+3918'
        cor1 = obs.Source.get_coordinates_from_name(src_name1)
        cor2 = obs.Source.get_coordinates_from_name(src_name2)
        assert cor1 == cor2
        assert cor1.to_string('hmsdms') == '00h00m41.527583s +39d18m04.14836s'
        src1 = obs.Source.get_source_from_name(src_name1)
        src2 = obs.Source.get_source_from_name(src_name2)
        assert src1.coord == src2.coord
        assert src1.coord == cor2
    except coord.name_resolve.NameResolveError:
        rprint("[yellow]WARNING: seems like no internet connection is available. " \
               "Unable to resolve source names[/yellow]")



def test_polarizations():
    assert obs.Polarization.SINGLE == obs.Polarization(1)
    assert obs.Polarization.SINGLE == obs.Polarization['SINGLE']
    assert obs.Polarization.DUAL == obs.Polarization(2)
    assert obs.Polarization.DUAL == obs.Polarization['DUAL']
    assert obs.Polarization.FULL == obs.Polarization(4)
    assert obs.Polarization.FULL == obs.Polarization['FULL']


def test_observation_init():
    o = obs.Observation()
    o.polarizations = 'dual'
    assert o.polarizations.value == 2
    o.polarizations = 4
    assert o.polarizations == obs.Polarization.FULL
    target = obs.Source('Target', '10h58m29.6s +81d33m58.8s')

    o = obs.Observation(target=target)
    o.times = Time('2020-06-15 20:00', scale='utc') + np.arange(0, 720, 10)*u.min
    o.band = '18'
    o.datarate = 1024
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

    o.inttime = 2

    # all_stations = stations.Network.get_stations_from_configfile()
    evn6 = ['Ef', 'Jb2', 'On', 'Hh', 'T6', 'Wb', 'Sv', 'Zc', 'Pa', 'Mp', 'Ho', 'Nl', 'Pt', 'Sc', 'Kp', 'Hn']
    o.stations_from_codenames(evn6)
    _ = o.elevations()
    _ = o.is_visible()
    # beam = o.synthesized_beam()
    # _ = o.thermal_noise()
    # _ = o.get_uv_array()
    # _ = o.get_dirtymap(robust='natural')
    # _ = o.get_dirtymap(robust='uniform')
    # fig, ax = plt.subplots()


