import pytest
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
