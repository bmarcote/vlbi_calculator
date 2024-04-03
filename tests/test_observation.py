import pytest
import numpy as np
from vlbiplanobs import observation as obs
from astropy import coordinates as coord
from astropy import units as u
from astropy.time import Time


def test_source():
    """Tests the Source class
    """
    with pytest.raises(AssertionError):
        obs.Source()

    with pytest.raises(ValueError):
        obs.Source('coordinates', 'a_name')
        obs.Source('00h00m00s xxhqqmees', 'a_name')
        obs.Source('40h00m00s 10d10m10s', 'a_name')
        obs.Source('00:00:00 100:00:00', 'a_name')
        obs.Source('-30:00:00 10:00:00', 'a_name')
        obs.Source('WQQWE')
        obs.Source('Cyg X-1')

    with pytest.raises(coord.name_resolve.NameResolveError):
        obs.Source(name='WQQWE')

    s1 = obs.Source('10h20m10s 40d30m10s', 'a_name')
    s2 = obs.Source('10:20:10 40:30:10', 'a_name', unit=(u.hourangle, u.deg))
    s3 = obs.Source(name='Cyg X-1')
    assert s1.name == 'a_name'
    assert s2.name == 'a_name'
    assert s1.coord == s2.coord
    assert s3.name == 'Cyg X-1'

    # checking that the distance to the Sun is right
    # times are spring equinox, summer solstice fall equinox
    times = Time(['2024-03-19 22:06', '2024-06-20 15:51', '2024-09-22 07:44'], scale='utc')
    s1 = obs.Source('0h0m0s 0d0m0s', 'equinox source')
    # print(f"\n\nSeparations: {s1.sun_separation(times)}\n\n")
    assert all(np.abs(s1.sun_separation(times) - np.array([0.0, 90.0, 180])*u.deg) < 1*u.deg)
    # print(f"\n\nConstraints: {s1.sun_constraint(20*u.deg, times=times)}\n\n")
    assert len(s1.sun_constraint(20*u.deg, times=times)) == 1
    s1 = obs.Source('0h0m0s 90d0m0s', 'polar source')
    print(f"\n\nSeparations: {s1.sun_separation(times)}\n\n")
    assert all(np.abs(s1.sun_separation(times[(0, 2),]) - np.array([90.0, 90.0])*u.deg) < 1*u.deg)
