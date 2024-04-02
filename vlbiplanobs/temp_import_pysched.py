# -*- coding: utf-8 -*-
# Licensed under GPLv3+ - see LICENSE

import configparser
from importlib import resources





def read_station_file():
    config = configparser.ConfigParser()
    with resources.as_file(resources.files("data").joinpath("stations_catalog.inp")) \
                                                              as stations_catalog_path:
        config.read(stations_catalog_path)

    return config
    # if all([key in config[stationname] for key in ('mount', 'ax1rate', 'ax2rate', 'ax1lim',
    #                                                    'ax2lim')]):


def read_pysched_stations(path='/Users/hawky/.pysched/catalogs/stations_RDBE.dat'):
    with open(path, 'r') as sched:
        all_lines = []
        for aline in sched.readlines():
            if (len(aline.strip()) > 0) and (aline.strip()[0] != '!'):
                all_lines.append(aline)

        entries = [e.strip() for e in ' '.join(all_lines).replace('\n', ' ').split('/')]
        # print(f"First entry: {entries[0]}")
        good_entries = []
        for i,entry in enumerate(entries):
            station_entry = []
            temp = [e.strip() for e in entry.split() if e.strip() != '']
            # print(f"Temp entry: {temp}")
            for t in temp:
                # print(f"T entry: {t}")
                if ('=' in t) and ('=' != t):
                    station_entry.append(t)
                elif len(station_entry) == 0:
                    station_entry.append(t)
                elif ('=' in station_entry[-1]) and ('=' != station_entry[-1][-1]):
                    station_entry.append(t)
                else:
                    station_entry[-1] += t

            good_entries.append(station_entry)

        all_antennas = {}
        # print(f"Good entries: {good_entries[0:4]}")
        for entry in good_entries:
            # print(f"\n\nENTRY: {entry}\n\n")
            ant = {e.split('=')[0].strip(): e.split('=')[1].strip() for e in entry if '=' in e}
            # print(ant)
            if ant != {}:
                all_antennas[ant['STCODE']] = ant

        return all_antennas



def add_sched2_stations():
    mystations = read_station_file()
    schedstats = read_pysched_stations()

    # print(f"Sched stations: {', '.join(schedstats)}")
    for ant in mystations.sections():
        if 'mount' not in mystations[ant]:
            if mystations[ant]['code'] in schedstats:
                for key in ('mount', 'ax1lim', 'ax2lim', 'ax1rate', 'ax2rate', 'ax1acc', 'ax2acc'):
                    if key.upper() in schedstats[mystations[ant]['code']]:
                        mystations[ant][key] = schedstats[mystations[ant]['code']][key.upper()]
            else:
                print(f"WARNING: station {mystations[ant]['code']} not found in pySCHED catalog.")

    with resources.as_file(resources.files("data").joinpath("stations_catalog.inp")) \
                                                                                as newfile:
        with open(newfile, 'w') as outfile:
            config = configparser.ConfigParser()
            config = mystations
            # config.read(newfile)
            config.write(outfile)



if __name__ == '__main__':
    add_sched2_stations()



