import sys
import argparse
from typing import Optional
import numpy as np
from astropy import units as u
from astropy.time import Time
from rich import print as rprint
from rich_argparse import RawTextRichHelpFormatter
import plotext as pltt
from vlbiplanobs import stations as stats
from vlbiplanobs import observation as obs
from vlbiplanobs import sources as src
# from vlbiplanobs import scheduler

_STATIONS = stats.Stations()
_NETWORKS = _STATIONS.get_networks_from_configfile()


def get_stations(list_networks: Optional[list[str]] = None,
                 list_stations: Optional[list[str]] = None) -> stats.Stations:
    """Returns a VLBI array including the required stations.
    Each argument is a comma-separated list of names.

    Inputs
        list_networks : list[str] | None
            If you want to pick the default antennas participating in one of the
            known VLBI networks, then you can include directly the network name
            here.
        list_stations : list[str] | None
            If you want a particular list of stations, or adding some that are not
            in the default network, then you can quote them here, using either the
            station code names or their names.
    """
    stations = []
    try:
        if list_networks is not None:
            networks = [_NETWORKS[n] for n in list_networks]
            for n in networks:
                for s in n.station_codenames:
                    if s not in stations:
                        stations.append(s)
    except KeyError:
        unknown_networks: list = [n for n in list_networks if n not in _NETWORKS]
        n = len(unknown_networks)
        rprint(f"[bold red]The network{'s' if n > 1 else ''} {', '.join(unknown_networks)}"
               f" {'are' if n > 1 else 'is'} not known.[/bold red]")
        sys.exit(1)

    try:
        if list_stations is not None:
            for s in list_stations:
                a_station = _STATIONS[s.strip()].codename
                if a_station not in stations:
                    stations.append(a_station)
    except KeyError:
        rprint(f"[bold red]The station {a_station} is not known.[/bold red]")
        sys.exit(1)

    return _STATIONS.filter_antennas(stations)


def summary(o: obs.Observation):
    """Prints the infromation for the given Observation, for testing purposes
    """
    rprint("\n[bold green]VLBI observation[/bold green]")
    rprint(f"To be conducted at {o.band.replace('cm', ' cm')} ", end='')
    if o.times is not None:
        rprint(f"from {o.times[0].strftime('%d %b %Y %H:%M')}--"
               f"{o.times[-1].strftime('%H:%M')} UTC")
    else:
        rprint("at unspecified times.")

    if o.duration is not None:
        rprint(f"With a total duration of {o.duration.to(u.h)}.")

    rprint(f"\n[bold green]Setup[/bold green]")
    if None not in (o.datarate, o.bandwidth, o.subbands):
        rprint(f"\nData rate of {o.datarate}, "
               f"{o.subbands} x {int(o.bandwidth.value/o.subbands)} "
               f"{o.bandwidth.unit} subbands, with {o.channels} channels each, "
               f"{o.polarizations} polarization.")
    else:
        rprint("[dim]No setup (data rate, bandwidth, number of subbands) specified[/dim]")

    rprint(f"\n[bold green]Stations ({len(o.stations)})[/bold green]: "
           f"{', '.join(o.stations.station_codenames)}")
    rprint("\n[bold green]Sources[/bold green]:")
    if o.scans is not None:
        for ablockname, ablock in o.scans.items():
            rprint(f"    - [dim]ScanBlock[/dim] '{ablockname}'\n      "
                   f"{'\n      '.join([s.name + ' [dim](' +
                      s.coord.to_string('hmsdms') + ')[/dim]'
                      for s in ablock.sources()])}")

    print('\n')


def plot_visibility(o: Observation, min_stations: int = 5):
    """Show plots with the different sources and when they are visible within the
    observation
    """
    elevs = o.elevations()
    srcup = o.is_observable()
    srcupalways = o.is_always_observable()
    when = o.when_is_observable(min_stations=5)

def main(band: str, networks: Optional[list[str]] = None,
         stations: Optional[list[str]] = None,
         src_catalog: Optional[str] = None, targets: Optional[list[str]] = None,
         start_time: Optional[Time] = None,
         duration: Optional[u.Quantity | float] = None):
    """Planner for VLBI observations.

    Inputs
        targets: list[str]
            List of sources to be observed. Each entry can be:
            a) the name defining a block in the source catalog file (if provided),
            b) the coordinates of the source, in RA, DEC (J2000), as 'hh:mm:ss dd:mm:ss'
               or 'XXhXXmXXs XXdXXmXXs',
            c) the name of the source, if it is a known one so it can be found in the
               SIMBAD/NEW/VizieR databases.
            A mix of the previous ones can also be used for each entry.
    """
    if networks is None and stations is None:
        rprint("[bold red]You need to provide at least a VLBI network "
               "or a list of antennas that will participate in the observation.[/bold red]")
        sys.exit(1)

    if isinstance(duration, float):
        duration = duration * u.hour

    if start_time is not None and start_time.scale != 'utc':
        rprint("[bold red]The start time must be in UTC[/bold red]\n"
               "[red](use 'scale' when defining the Time object)[/red]")
        sys.exit(1)

    if src_catalog is not None:
        source_catalog = src.SourceCatalog(src_catalog)

    src2observe: dict[str, src.ScanBlock] = {}
    if targets is not None:
        for target in targets:
            if target in source_catalog.blocknames:
                src2observe[target] = source_catalog[target]
            else:
                a_source = src.Source.source_from_str(target)
                src2observe[a_source.name] = src.ScanBlock([src.Scan(a_source)])
    elif src_catalog is None:
        rprint("[bold red]Either a source catalog file or a list of targets must be "
               "provided (or both).[/bold red]")
        sys.exit(1)
    else:
        src2observe = source_catalog.blocks

    o = obs.Observation(band, get_stations(networks, stations), scans=src2observe,
                        times=start_time + np.arange(0, duration.to(u.min).value, 10)*u.min
                        if start_time is not None else None, duration=duration,
                        # datarate=, subbands=, channels=, polarizations=, inttime=)
                        ontarget=0.6)
    summary(o)


if __name__ == '__main__':
    usage = "%(prog)s [-h]  XXX"
    description = "EVN Observation Planner"
    parser = argparse.ArgumentParser(description=description, prog="planobs", usage=usage,
                                     formatter_class=RawTextRichHelpFormatter)
    parser.add_argument('-t', '--targets', type=str, default=None,
                        help="Source(s) to be observed.")
    parser.add_argument('-i', '--input', type=str, default=None,
                        help="Input file containing the personal source catalog. "
                        "If provided, then the '--targets' will select the block(s) "
                        "defined in this file.")
    parser.add_argument('-t1', '--starttime', type=str, default=None,
                        help="Start of the observation, with the shape YYYY/MM/DD/HH:MM "
                        "in UTC.")
    parser.add_argument('-d', '--duration', type=str, default=None,
                        help="Total duration of the observation, in hours.")
    parser.add_argument('-n', '--network', type=str, help="Comma-separated list "
                        "of the VLBI network(s) that will participate in the observation. "
                        "It will take the default stations in each network. If 'stations' "
                        "is provided, then it will take both the default stations plus "
                        "the ones given in stations.")
    parser.add_argument('-s', '--stations', type=str, help="Comma-separated list "
                        "of the antennas that will participate in the observation. "
                        "You can use either antenna codenames or the standard name "
                        "as given in the catalogs. See '--list-stations' to get a list.")
    parser.add_argument('-b', '--band', type=str, help="Observing band, as defined"
                        "in the catalogs as 'XXcm', with 'XX' the wavelegnth in cm. "
                        "See '--list-bands' to get a list.")
    # parser.add_argument('', type=, default=, help='')
    parser.add_argument('--list-antennas', action="store_true", default=False,
                        help="Writes the list of all antennas defined in PlanObs.")
    parser.add_argument('--list-bands', action="store_true", default=False,
                        help="Writes the list of all observing bands defined in PlanObs.")
    # parser.add_argument('', type=, default=, help='')
    subparser = parser.add_subparsers(dest='subparser', help="")
    parser.add_argument('--sched', default=None, type=str,
                        help="Produces a (SCHED) .key schedule file for "
                        "the observation with the given name.")
    parser.add_argument('--fringefinders', default='2', type=str,
                        help="Defines the fringe finder sources to be scheduled "
                        "in the observation. It can be either a comma-separated "
                        "list of source names (as long as they appear in AstroGeo) "
                        "or a number, meaning how many scans should go on fringe "
                        "finders, and it will select the most suitable sources.")
    parser.add_argument('--polcal', action="store_true", default=False,
                        help="Requires polarization calibration for the observation")
    parser.add_argument('--pulsar', default=None, type=str,
                        help="Sets to schedule at least a scan on a pulsar source. "
                        "If a number, it will select a pulsar from the personal "
                        "input source file (must be provided!). If a name, "
                        "it will pick such pulsar.")

    args = parser.parse_args()

    if args.list_antennas:
        rprint("[bold]Available VLBI networks:[/bold]")
        for network_name, network in _NETWORKS.items():
            rprint(f"[bold]{network_name}[/bold]")
            rprint(f"  [dim]Default antennas: {', '.join(network.station_codenames)}[/dim]")

        rprint("\n[bold]All available antennas:[/bold]")
        for ant in _STATIONS:
            rprint(f"    {ant.codename} ({ant.name})  {ant.diameter} in {ant.country}")
            rprint(f"      [dim]Observes at {', '.join(ant.bands)}[/dim]")

    if args.list_bands:
        rprint("\n[bold]Available observing bands:[/bold]")
        for aband in obs.freqsetups.bands:
            rprint(f"[bold]{aband}[/bold] [dim]({obs.freqsetups.bands[aband]})[/dim]")
            rprint(f"[dim]  Observable with {', '.join([nn for nn, n in _NETWORKS.items()
                                             if aband in n.observing_bands])}[/dim]")

    if args.list_antennas or args.list_bands:
        sys.exit(0)

    if args.band not in obs.freqsetups.bands:
        rprint(f"[bold red]The provided band ({args.band}) is not available"
               "[/bold red]\n[dim]Available values are: "
               f"{', '.join(obs.freqsetups.bands)}.[/dim]")
        sys.exit(1)

    if args.network is None and args.stations is None:
        rprint("[bold red]You need to provide at least a VLBI network "
               "or a list of antennas that will participate in the observation.[/bold red]")
        sys.exit(1)

    main(band=args.band,
         networks=[n.strip() for n in args.network.split(',')] if args.network else None,
         stations=[s.strip() for s in args.stations.split(',')] if args.stations else None,
         src_catalog=args.input,
         targets=args.targets.split(','), start_time=Time(args.starttime, scale='utc')
         if args.starttime else None, duration=args.duration*u.hour if args.duration else None)
    # o = obs.Observation(band=args.band, stations=get_stations(args.network, args.stations))
