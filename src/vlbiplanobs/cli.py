









# This is code taken from the Observation class, likely to be integrated in the potential VLBIPlanObs class
# to be used in a command line interface.  And should merge stations and sources so it makes the user's life
# much easier.



    def with_networks(self, networks: Union[str, list[str]]):
        """Selects the given network(s) to observe in this observation. It can be only one single
        VLBI network or multiple of them.
        This function will adjust the setup of the observations to the default setup for the given networks,
        although this can be later manually adjust again by the user.

        Input
            networks :  str | list[str]
                Name fo the network or networks to join the observation.
        """
        if isinstance(networks, str):
            self.stations = _NETWORKS[networks]
        elif isinstance(networks, list):
            self.stations = _NETWORKS[networks[0]]
            self.stations.stations = [s for n in networks for s in _NETWORKS[n].stations]
        else:
            raise ValueError(f"The value of 'networks' needs to be a str or a list of str, not '{networks}'.")

    def stations_from_codenames(self, new_stations: list[str]):
        """Includes the given antennas in the observation, specified by their code names as
        set in the stations_catalog file.
        """
        self.stations = Network.get_stations_from_configfile(codenames=new_stations)

    # TODO: re-run this function if the band changes!
    def sources_from_catalog(self, path: str):
        """Reads the yaml file with the information of all sources that may be scheduled.
        It returns the catalog as a dictionary.

        Input
            path : str
                The path to the yaml file with the catalog of sources to be imported.
        """
        with open(path, 'r') as sources_yaml:
            catalog = yaml.safe_load(sources_yaml)
            scanblocks = []
            for a_type in catalog:
                for a_src in catalog[a_type]:
                    scans = []
                    src = catalog[a_type][a_src]
                    if 'phasecal' in src:
                        scans.append(Scan(Source(name=src['phasecal']['name'],
                                                 coordinates=src['phasecal']['coordinates'],
                                                 source_type=SourceType.PHASECAL),
                                          duration=src['phasecal']['duration']*u.min
                                          if 'duration' in src['phasecal'] else 2*u.min))

                    if 'checkSource' in src:
                        scans.append(Scan(Source(name=src['checkSource']['name'],
                                                 coordinates=src['checkSource']['coordinates'],
                                                 source_type=SourceType.CHECKSOURCE),
                                          duration=src['checkSource']['duration']*u.min
                                          if 'duration' in src['checkSource'] else 3*u.min,
                                          every=src['checkSource']['repeat']
                                          if 'repeat' in src['checkSource'] else 3))

                    match a_type:
                        case 'pulsars':
                            src_type = SourceType.PULSAR
                        case 'targets':
                            src_type = SourceType.TARGET
                        case 'amplitudecals':
                            src_type = SourceType.AMPLITUDECAL
                        case 'checksources':
                            src_type = SourceType.CHECKSOURCE
                        case 'phasecals':
                            src_type = SourceType.PHASECAL
                        case 'fringefinders':
                            src_type = SourceType.FRINGEFINDER
                        case 'polcals':
                            src_type = SourceType.POLCAL
                        case _:
                            src_type = SourceType.UNKNOWN

                    if 'duration' not in src:
                        if (recomm := self.recommended_phaseref_cycle()) is not None:
                            src['duration'] = max(recomm - 2*u.min, 1*u.min)
                        else:
                            src['duration'] = 3.5*u.min
                    else:
                        src['duration'] *= u.min

                    scans.append(Scan(source=Source(name=a_src, coordinates=src['coordinates'],
                                                    source_type=src_type), duration=src['duration']))

                    scanblocks.append(ScanBlock(scans))
        self.scans = scanblocks
        # return scanblocks





