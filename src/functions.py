"""Different functions that are required to operate the program
"""

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
        if a_line['station'] is 'Arecibo':
            min_elev = 80*u.deg
        else:
            min_elev = 10*u.deg

        new_station = stations.SelectedStation(a_line['station'], a_line['code'],
             network=a_line['network'], location=a_loc, freqs_sefds=sefds,
             min_elevation=min_elev, selected=False)

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


