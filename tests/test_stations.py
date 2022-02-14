import pytest
import numpy as np
from vlbiplanobs import stations
from astropy import coordinates as coord
from astropy import units as u
from astropy.time import Time
from astroplan import FixedTarget


def test_station_init():
    """Tests the Station class and different possibilities during the creation of the station.
    """
    sefds = {'18': 100, '6': 40, '0.1': 200}
    a_station = stations.Station('name', 'Nm', 'VLBI', coord.EarthLocation(0., 0., 0.), sefds, 10)
    a_station = stations.Station('name', 'Nm', 'VLBI', coord.EarthLocation(0., 0., 0.), {}, 10.4)
    assert isinstance(a_station.min_elevation, u.Quantity)
    a_station = stations.Station('name', 'Nm', 'VLBI', coord.EarthLocation(0., 0., 0.), {}, 10.3*u.deg)
    with pytest.raises(AssertionError):
        # Wrong names
        a_station = stations.Station(None, 'Nm', 'VLBI', coord.EarthLocation(0., 0., 0.), {}, 10)
        a_station = stations.Station(5, 'Nm', 'VLBI', coord.EarthLocation(0., 0., 0.), {}, 10)
        a_station = stations.Station('name', None, 'VLBI', coord.EarthLocation(0., 0., 0.), {}, 10)
        a_station = stations.Station('name', True, 'VLBI', coord.EarthLocation(0., 0., 0.), {}, 10)
        a_station = stations.Station('name', 'Nm', None, coord.EarthLocation(0., 0., 0.), {}, 10)
        a_station = stations.Station('name', 'Nm', 5.0, coord.EarthLocation(0., 0., 0.), {}, 10)
        # Negative min_elevations
        a_station = stations.Station('name', 'Nm', 'VLBI', coord.EarthLocation(0., 0., 0.), {}, -10)
        # Wrong units min_elevations
        a_station = stations.Station('name', 'Nm', 'VLBI', coord.EarthLocation(0., 0., 0.), {}, 10*u.m)


def test_station_functions():
    """Tests the Station functions.
    """
    sefds = {'18': 100, '6': 40, '0.1': 200}
    a_station = stations.Station('name', 'Nm', 'VLBI',
                                 coord.EarthLocation(3839348.973*u.m, 430403.51*u.m, 5057990.099*u.m), sefds, 20)
    assert isinstance(a_station.name, str)
    assert isinstance(a_station.fullname, str)
    assert a_station.fullname == a_station.name
    assert a_station.network == 'VLBI'
    assert a_station.all_networks == a_station.network
    assert isinstance(a_station.country, str)
    assert isinstance(a_station.diameter, str)
    assert isinstance(a_station.real_time, bool)
    assert isinstance(a_station.location, coord.EarthLocation)
    assert a_station.location == coord.EarthLocation(3839348.973*u.m, 430403.51*u.m, 5057990.099*u.m)
    assert list(a_station.bands) == ['18', '6', '0.1']
    assert isinstance(a_station.sefds, dict)
    assert a_station.sefds == sefds
    assert a_station.has_band('18')
    assert not a_station.has_band('45')
    assert a_station.sefd('18') == 100
    with pytest.raises(KeyError):
        a_station.sefd('45')

    times1 = Time('2020-03-21 1:00') + np.arange(0, 4*60, 10)*u.min
    times2 = Time('2020-03-21 3:00') + np.arange(0, 4*60, 10)*u.min
    src1 = FixedTarget(coord=coord.SkyCoord('0h0m0s 30d0m0s'), name='testSrc')
    # a_station.elevation(times1, src1)  # Should have elevation ranging -5 to 16.8 deg.
    # a_station.elevation(times1, src1)  # Should have elevation ranging 3.7 to 34 deg.
    assert len(a_station.is_visible(times1, src1)[0]) == 0
    assert len(a_station.is_visible(times2, src1)[0]) == 10
    assert len(a_station.elevation(times2, src1)) == len(times1)
    assert np.equal(a_station.elevation(times2, src1).value, a_station.altaz(times2, src1).alt.value)[0]


