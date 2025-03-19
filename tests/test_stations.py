import pytest
import configparser
import numpy as np
from importlib import resources
from astropy import coordinates as coord
from astropy import units as u
from astropy.time import Time
from astroplan import FixedTarget
from vlbiplanobs import stations
from vlbiplanobs import freqsetups


def test_stations_catalog():
    """Checks that the stations_catalog.inp file is self-consistent.
    So no duplicated stations (with the same stationname, as that must be unique), or wrong parameters.
    """
    config = configparser.ConfigParser()
    with resources.as_file(resources.files("vlbiplanobs.data").joinpath("stations_catalog.inp")) as stat_cat:
        config.read(stat_cat)

    codenames = set()
    for stationname in config.sections():
        for a_key in ('station', 'code', 'networks', 'country', 'diameter',
                      'real_time', 'position'):
            assert a_key in config[stationname], \
                   f"'{a_key}' is not present in the stations_catalog.inp for {stationname}"

        assert config[stationname]['code'] not in codenames, \
               f"Duplicated codename {config[stationname]['code']} in stations_catalog.inp."
        codenames.add(config[stationname]['code'])
        temp = [float(i.strip()) for i in config[stationname]['position'].split(',')]
        assert len(temp) == 3, f"The station position of {stationname} does not have the three expected numbers."
        _ = coord.EarthLocation(u.Quantity(float(temp[0]), u.m), u.Quantity(float(temp[1]), u.m),
                                u.Quantity(float(temp[2]), u.m))
        # Getting the SEFD values for the bands
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
    with resources.as_file(resources.files("vlbiplanobs.data").joinpath("network_catalog.inp")) as net_cat_path:
        config.read(net_cat_path)

    # To verify that all stations (codenames) defined in the network file as defined in the stations_catalog
    all_stations = stations.Stations().station_codenames
    for networkname in config.sections():
        for a_key in ('name', 'default_antennas', 'max_datarate', 'observing_bands'):
            assert a_key in config[networkname], \
                   f"'{a_key}' is not present in the network_catalog.inp for {networkname}"
        network_antennas = [a.strip() for a in config[networkname]['default_antennas'].split(',')]
        for antenna in network_antennas:
            assert antenna in all_stations, \
                   f"{antenna}, defined in {networkname} is not present in stations_catalog.inp"


def test_station_init():
    """Tests the Station class and different possibilities during the creation of the station.
    """
    sefds = {'18cm': 100, '6cm': 40, '0.1cm': 200}
    _ = stations.Station('name', 'Nm', 'VLBI', coord.EarthLocation(0., 0., 0.), sefds)
    _ = stations.Station('name', 'Nm', ('EVN', 'VLBA'), coord.EarthLocation(0., 0., 0.), {})
    with pytest.raises(AssertionError):
        # Wrong names
        _ = stations.Station(None, 'Nm', 'VLBI', coord.EarthLocation(0., 0., 0.), {}, 10)
        _ = stations.Station(5, 'Nm', 'VLBI', coord.EarthLocation(0., 0., 0.), {}, 10)
        _ = stations.Station('name', None, 'VLBI', coord.EarthLocation(0., 0., 0.), {}, 10)
        _ = stations.Station('name', True, 'VLBI', coord.EarthLocation(0., 0., 0.), {}, 10)
        _ = stations.Station('name', 'Nm', None, coord.EarthLocation(0., 0., 0.), {}, 10)
        _ = stations.Station('name', 'Nm', 5.0, coord.EarthLocation(0., 0., 0.), {}, 10)
        # Negative min_elevations
        _ = stations.Station('name', 'Nm', 'VLBI', coord.EarthLocation(0., 0., 0.), {}, -10)
        # Wrong units min_elevations
        _ = stations.Station('name', 'Nm', 'VLBI', coord.EarthLocation(0., 0., 0.), {}, 10*u.m)


# TODO: add tests on: check Mount type is what is written in the file, and depending on that, if it got the
# correct constraints (e.g. HourAngle for EQUAT, etc)

def test_station_functions():
    """Tests the Station functions.
    """
    sefds = {'18cm': 100*u.Jy, '6cm': 40*u.Jy, '0.1cm': 200*u.Jy}
    a_station = stations.Station('name', 'Nm', 'VLBI',
                                 coord.EarthLocation(3839348.973*u.m, 430403.51*u.m, 5057990.099*u.m), sefds)
    assert isinstance(a_station.name, str)
    assert isinstance(a_station.fullname, str)
    assert a_station.fullname == a_station.name
    assert a_station.networks == ('VLBI',)
    assert len(a_station.networks) == 1
    assert isinstance(a_station.country, str)
    assert isinstance(a_station.diameter, str)
    assert isinstance(a_station.real_time, bool)
    assert isinstance(a_station.location, coord.EarthLocation)
    assert a_station.location == coord.EarthLocation(3839348.973*u.m, 430403.51*u.m, 5057990.099*u.m)
    assert list(a_station.bands) == ['18cm', '6cm', '0.1cm']
    assert isinstance(a_station.sefds, dict)
    assert a_station.sefds == sefds
    assert a_station.has_band('18cm')
    assert not a_station.has_band('45cm')
    assert a_station.sefd('18cm') == 100*u.Jy
    with pytest.raises(KeyError):
        a_station.sefd('45')

    times1 = Time('2020-03-21 1:00') + np.arange(0, 2*60, 10)*u.min
    times2 = Time('2020-03-21 3:00') + np.arange(0, 4*60, 10)*u.min
    src1 = FixedTarget(coord=coord.SkyCoord('0h0m0s 30d0m0s'), name='testSrc')
    # a_station.elevation(times1, src1)  # Should have elevation ranging -7 to 1.1 deg.
    # a_station.elevation(times2, src1)  # Should have elevation ranging 2.2 to 34 deg.
    assert not all(a_station.is_observable(times1, src1))
    assert all(a_station.is_observable(times2, src1))
    assert not a_station.is_always_observable(times2, src1)
    assert len(a_station.elevation(times2, src1)) == len(times2)
    assert len(a_station.elevation(times1, src1)) == len(times1)
    assert np.equal(a_station.elevation(times2, src1).value, a_station.altaz(times2, src1).alt.value)[0]


def test_station_file():
    """Tests the consistency between the antennas in the station catalog, network catalog, and freq setups.
    Mainly:  no duplicated antenna code names, all antennas in the networks need to be defined, and all
    setups in the antennas need to be defined in the freqsetup.py
    """
    # all_stations = list(stations.Stations.get_stations_from_configfile())
    all_stations = stations.Stations()
    assert len(all_stations) > 0

    # Verifying networks - it would raise issues if inconsistenties
    all_networks = stations.Stations.get_networks_from_configfile()

    # double checking with one Station in particular
    with open(resources.files("vlbiplanobs.data").joinpath("stations_catalog.inp")) as afile:
        in_ef = False
        n_ants = 0
        for aline in afile.readlines():
            if 'station=' in aline.strip().replace(' ', ''):
                n_ants += 1

            if in_ef and ('station=' in aline.strip().replace(' ', '')):
                # We finished with effelsberg!
                in_ef = False

            if in_ef and (len(aline.strip()) > 0) and (aline.strip()[0] != '#') and ('[' not in aline):
                key, value = [a.strip() for a in aline.split('=')]
                match key:
                    case 'station':
                        assert value == all_stations['Ef'].fullname
                    case 'network':
                        assert value == all_stations['Ef'].network
                    case 'possible_networks':
                        for avalue in value.split(','):
                            assert avalue.strip() in all_stations['Ef'].all_networks
                    case 'country':
                        assert value == all_stations['Ef'].country
                    case 'diameter':
                        assert value == all_stations['Ef'].diameter
                    case 'position':
                        avalue = [float(a.strip()) for a in value.split(',')]
                        assert avalue[0]*u.m == all_stations['Ef'].location.x
                        assert avalue[1]*u.m == all_stations['Ef'].location.y
                        assert avalue[2]*u.m == all_stations['Ef'].location.z
                    case 'mount':
                        assert value == all_stations['Ef'].mount.mount_type.name
                    case 'ax1rate':
                        assert float(value)*u.deg/u.s == all_stations['Ef'].mount.ax1.speed
                    case 'ax2lim':
                        avalue = [float(a.strip())*u.deg for a in value.split(',')]
                        assert avalue[0] == all_stations['Ef'].mount.ax2.limits[0]
                        assert avalue[1] == all_stations['Ef'].mount.ax2.limits[1]
                    case 'ax1acc':
                        assert float(value)*u.deg/u.s/u.s == all_stations['Ef'].mount.ax1.acceleration
                    case 'real_time':
                        assert bool(value) == all_stations['Ef'].real_time

                if 'SEFD' in key:
                    assert all_stations['Ef'].has_band(key.removeprefix('SEFD_') + 'cm')
                    assert all_stations['Ef'].sefd(key.removeprefix('SEFD_')) == float(value)
                    assert key.removeprefix('SEFD_') + 'cm' in freqsetups.bands

            if aline.strip() == 'station = Effelsberg':
                in_ef = True

            # ingore if we are at any other station, but:
            # At the same time, I check for consistency in the network names
            if aline.strip().replace(' ', '') == 'network=':
                value = aline.split('=')[1].strip()
                assert ',' not in value
                assert value in all_networks, "The network '{value}' is not defined in the networks file."

            if aline.strip().replace(' ', '') == 'possible_network=':
                values = [a.strip() for a in aline.split('=')[1].strip().split(',')]
                for avalue in values:
                    assert avalue in all_networks, "The network '{avalue}' is not defined in the networks file."

            if aline.strip().replace(' ', '') == 'maxdatarate_':
                assert aline.split('=')[0].strip().removeprefix('maxdatarate_') in all_networks

    # This value may change in the future but then is a problem of the test, not the code
    assert 'EVN' in all_stations['De'].max_datarate, \
           f"De doesn't have max_datarate for EVN (keys: {all_stations['De'].max_datarate})"
    assert all_stations['De'].max_datarate['EVN'] == 512*u.Mb/u.s, \
           f"De has a mismatching data rate of {all_stations['De'].max_datarate['EVN']}."
    assert n_ants == len(all_stations)

    # Check that all stations have an associated image
    for ant in all_stations:
        with resources.as_file(resources.files("vlbiplanobs.gui.assets").
                               joinpath(f"ant-{ant.name.replace(' ', '_').lower()}.jpg")) as antfile:
            assert antfile.exists()

    for net in all_networks:
        with resources.as_file(resources.files("vlbiplanobs.gui.assets").
                               joinpath(f"network-{net.replace(' ', '_').lower()}.png")) as netfile:
            if net != 'e-EVN':
                assert netfile.exists()
