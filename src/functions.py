"""Different functions that are required to operate the program
"""

import configparser
from astropy import coordinates as coord
from astropy import units as u
from astropy.io import ascii

from src import stations


def get_networks_from_file(filename='data/station_location.txt'):
    """Retrieves all the stations given in the filename and
    creates different arrays. Returns a dict with all arrays.
    The file must contain the information for each antenna on each
    row, with the following header (written at the beginning of the file):
        - station : full name of the station
        - code : code name of the station
        - network : name of the network the station belongs to.
        - x : x coordinates of the station, in meters.
        - y : y coordinates of the station, in meters.
        - z : z coordinates of the station, in meters.
        - SEFD-XX : SEFD of the station at XX cm (-1 if this band cannot be observed)
    """
    datafile = ascii.read(filename)
    # Creates a `Stations` object for each network, for now with no stations
    network_labels = set(datafile['network'])
    networks = {}
    for network in network_labels:
        networks[network] = stations.Stations(network, [])

    sefd_names = datafile.colnames[6:]
    for a_line in datafile:
        # Create each station and add it to the correct noetwork
        sefds = {}
        for sefd_name in sefd_names:
            if a_line[sefd_name] > 0.0:
                sefds[f"{sefd_name.replace('SEFD-','').strip()}cm"] = a_line[sefd_name]

        a_loc = coord.EarthLocation(a_line['x']*u.m, a_line['y']*u.m, a_line['z']*u.m)
        # For now this is hard-coded
        if a_line['station'] is 'Arecibo':
            min_elev = 80*u.deg
        else:
            min_elev = 10*u.deg

        new_station = stations.SelectedStation(a_line['station'], a_line['code'],
             location=a_loc, freqs_sefds=sefds, min_elevation=min_elev, selected=False)

        networks[a_line['network']].add(new_station)

    return networks

def get_stations_from_file(filename='data/station_location.txt'):
    """Retrieves all the stations given in the filename and
    creates a Stations object containing all stations, which is returned.
    The file must contain the information for each antenna on each
    row, with the following header (written at the beginning of the file):
        - station : full name of the station
        - code : code name of the station
        - network : name of the network the station belongs to.
        - x : x coordinates of the station, in meters.
        - y : y coordinates of the station, in meters.
        - z : z coordinates of the station, in meters.
        - SEFD-XX : SEFD of the station at XX cm (-1 if this band cannot be observed)
    """
    datafile = ascii.read(filename)
    networks = stations.Stations('network', [])
    sefd_names = datafile.colnames[6:]
    for a_line in datafile:
        # Create each station and add it to the correct noetwork
        sefds = {}
        for sefd_name in sefd_names:
            if a_line[sefd_name] > 0.0:
                sefds[f"{sefd_name.replace('SEFD-','').strip()}cm"] = a_line[sefd_name]

        a_loc = coord.EarthLocation(a_line['x']*u.m, a_line['y']*u.m, a_line['z']*u.m)
        # For now this is hard-coded
        if a_line['code'] == 'Ar':
            min_elev = 65*u.deg
        else:
            min_elev = 10*u.deg

        new_station = stations.SelectedStation(a_line['station'], a_line['code'],
             network=a_line['network'], location=a_loc, freqs_sefds=sefds,
             min_elevation=min_elev, selected=False)

        networks.add(new_station)

    return networks


def get_stations_from_configfile(filename='data/stations_catalog.inp'):
    """Retrieves the information concerning all stations available in the 'filename'
    file. Creates a Stations object containing the stations and the information on it.
    The file must have a format readable by the Python ConfigParser.
    Each section will be named with the name of the station, and then it must have
    the following keys:
    station - full name of the station.
    code - codename for the station (typically two letters).
    network - main network to which it belongs to.
    possible_networks - all networks the station can participate in (including 'network')
    country - country where the station is located.
    diameter - string with the diameter of the station.
    position = x, y, z (in meters). Geoposition of the station.
    min_elevation (in degrees) - minimum elevation the station can observe.
    SEFD_**  - SEFD of the station at the **cm band. If a given band is not present,
                it is assumed that the station cannot observe it.
    [optional]
    img - a path to an image of the station.
    link - a url linking to the station page/related information.
    """
    config = configparser.ConfigParser()
    config.read(filename)
    networks = stations.Stations('network', [])
    for stationname in config.sections():
        temp = [float(i.strip()) for i in config[stationname]['position'].split(',')]
        a_loc = coord.EarthLocation(temp[0]*u.m, temp[1]*u.m, temp[2]*u.m)
        # Getting the SEFD values for the bands
        min_elev = float(config[stationname]['min_elevation'])*u.deg
        sefds = {}
        for akey in config[stationname].keys():
            if 'SEFD_' in akey.upper():
                sefds[f"{akey.upper().replace('SEFD_', '').strip()}cm"] = \
                                    float(config[stationname][akey])
        new_station = stations.SelectedStation(stationname, config[stationname]['code'],
                config[stationname]['network'], a_loc, sefds, min_elev,
                config[stationname]['station'], config[stationname]['possible_networks'],
                config[stationname]['country'], config[stationname]['diameter'])
        networks.add(new_station)

    return networks


def stations_with_band(networks, band, output_network_name=None):
    """For a given collection of networks or Stations, returns a Stations object
    including all stations that can observe at the given band.

        - networks : dict [name_network]: Stations or Stations
        - band : str
    """
    if output_network_name is None:
        output_network_name = f"Stations@{band}"

    antennas = stations.Stations(output_network_name, [])
    if isinstance(networks, dict):
        for network in networks:
            for station in networks[network]:
                if band in station.bands:
                    antennas.add(station)
    elif isinstance(networks, stations.Stations):
        for station in networks:
            if band in station.bands:
                antennas.add(station)
    else:
        raise ValueError(f"{networks} expected to be either dict of Stations type.")
    return antennas


def print_obs_times(obs, date_format='%d %b %Y'):
    """Given an observation, it returns the time range (starttime-endtime) in a smart
    way. If the observation lasts for less than one day it omits the end date:
            20 Jan 1971 10:00-20:00UT
    It also adds the GST range after that.

    Input:
        - obs : observation.Observation
            It must already have set the .times part with an array of astropy.Time times.
        - date_format : str [optional]
            Format for the date part (only the date part) of the string to represent
            the time range.
    Output:
        - printed_time : str
            A string showing the time-range of the observation.

    """
    gsttext = "{:02n}:{:02.2n}-{:02n}:{:02.2n}".format((obs.gstimes[0].hour*60) // 60,
                                              (obs.gstimes[0].hour*60) % 60,
                                              (obs.gstimes[-1].hour*60) // 60,
                                              (obs.gstimes[0].hour*60) % 60)
    if obs.times[0].datetime.date() == obs.times[-1].datetime.date():
        return "{}\n{}-{} UTC\nGST: {}".format(obs.times[0].datetime.strftime(date_format),
                                    obs.times[0].datetime.strftime('%H:%M'),
                                    obs.times[-1].datetime.strftime('%H:%M'), gsttext)
    elif (obs.times[-1] - obs.times[0]) < 24*u.h:
        return "{}\n{}-{} UTC (+1d)\nGST: {}".format(
                                    obs.times[0].datetime.strftime(date_format),
                                    obs.times[0].datetime.strftime('%H:%M'),
                                    obs.times[-1].datetime.strftime('%H:%M'), gsttext)
    else:
        return "{} {} to {} {} UTC\nGST: {}".format(
                                    obs.times[0].datetime.strftime(date_format),
                                    obs.times[0].datetime.strftime('%H:%M'),
                                    obs.times[-1].datetime.strftime(date_format),
                                    obs.times[-1].datetime.strftime('%H:%M'), gsttext)













