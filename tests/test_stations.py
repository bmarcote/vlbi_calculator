import pytest
import configparser
from importlib import resources
import numpy as np
from vlbiplanobs import stations
from astropy import coordinates as coord
from astropy import units as u
from astropy.time import Time
from astroplan import Observer, FixedTarget




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


def test_stations_catalog():
    """Checks that the stations_catalog.inp file is self-consistent.
    So no duplicated stations (with the same stationname, as that must be unique), or wrong parameters.
    """
    config = configparser.ConfigParser()
    with resources.path("data", "stations_catalog.inp") as stations_catalog_path:
        config.read(stations_catalog_path)

    codenames = set()
    for stationname in config.sections():
        for a_key in ('station', 'code', 'network', 'possible_networks', 'country', 'diameter', \
                      'min_elevation', 'real_time', 'position'):
            assert a_key in config[stationname], \
                   f"'{a_key}' is not present in the stations_catalog.inp for {stationname}"

        assert config[stationname]['code'] not in codenames, \
               f"Duplicated codename {config[stationname]['code']} in stations_catalog.inp."
        codenames.add(config[stationname]['code'])
        temp = [float(i.strip()) for i in config[stationname]['position'].split(',')]
        assert len(temp) == 3, f"The station position of {stationname} does not have the three expected numbers."
        a_loc = coord.EarthLocation(float(temp[0])*u.m, float(temp[1])*u.m, float(temp[2])*u.m)
        # Getting the SEFD values for the bands
        assert float(config[stationname]['min_elevation']) >= 0.0,f"min_elevation for {stationname} must be >= 0.0"
        assert config[stationname]['real_time'] in ('yes', 'no'), \
               f"'real_time' must be 'yes' or 'no' for {stationname}"
        # At least one SEFD per station
        has_sefd = False
        for akey in config[stationname].keys():
            has_sefd = ('SEFD_' in akey.upper()) or has_sefd

        assert has_sefd, f"Station {stationname} must have at least a SEFD value for one frequency."


def test_network_catalog():
    """Checks the consistency of the information in the  network_catalog.inp file.
    """
    config = configparser.ConfigParser()
    with resources.path("data", "network_catalog.inp") as networks_catalog_path:
        config.read(networks_catalog_path)

    # To verify that all stations (codenames) defined in the network file as defined in the stations_catalog
    all_stations = stations.Stations.get_stations_from_configfile().codenames
    for networkname in config.sections():
        for a_key in ('name', 'default_antennas', 'max_datarate', 'observing_bands'):
            assert a_key in config[networkname], \
                   f"'{a_key}' is not present in the network_catalog.inp for {networkname}"
        network_antennas = [a.strip() for a in config[networkname]['default_antennas'].split(',')]
        for antenna in network_antennas:
            assert antenna in all_stations, \
                   f"{antenna}, defined in {networkname} is not present in stations_catalog.inp"
        temp = int(config[networkname]['max_datarate'])






