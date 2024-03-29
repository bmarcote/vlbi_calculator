import configparser
from importlib import resources
from vlbiplanobs import stations
from astropy import coordinates as coord
from astropy import units as u


def test_stations_catalog():
    """Checks that the stations_catalog.inp file is self-consistent.
    So no duplicated stations (with the same stationname, as that must be unique), or wrong parameters.
    """
    config = configparser.ConfigParser()
    with resources.as_file(resources.files("data").joinpath("stations_catalog.inp")) as stations_catalog_path:
        config.read(stations_catalog_path)

    codenames = set()
    for stationname in config.sections():
        for a_key in ('station', 'code', 'network', 'possible_networks', 'country', 'diameter', \
                      'real_time', 'position'):
            assert a_key in config[stationname], \
                   f"'{a_key}' is not present in the stations_catalog.inp for {stationname}"

        assert config[stationname]['code'] not in codenames, \
               f"Duplicated codename {config[stationname]['code']} in stations_catalog.inp."
        codenames.add(config[stationname]['code'])
        temp = [float(i.strip()) for i in config[stationname]['position'].split(',')]
        assert len(temp) == 3, f"The station position of {stationname} does not have the three expected numbers."
        a_loc = coord.EarthLocation(u.Quantity(float(temp[0]), u.m),        u.Quantity(float(temp[1]), u.m),u.Quantity(float(temp[2]), u.m))
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
    with resources.as_file(resources.files("data").joinpath("network_catalog.inp")) \
                                                                                         as networks_catalog_path:
        config.read(networks_catalog_path)

    # To verify that all stations (codenames) defined in the network file as defined in the stations_catalog
    all_stations = stations.Network.get_stations_from_configfile().codenames
    for networkname in config.sections():
        for a_key in ('name', 'default_antennas', 'max_datarate', 'observing_bands'):
            assert a_key in config[networkname], \
                   f"'{a_key}' is not present in the network_catalog.inp for {networkname}"
        network_antennas = [a.strip() for a in config[networkname]['default_antennas'].split(',')]
        for antenna in network_antennas:
            assert antenna in all_stations, \
                   f"{antenna}, defined in {networkname} is not present in stations_catalog.inp"


