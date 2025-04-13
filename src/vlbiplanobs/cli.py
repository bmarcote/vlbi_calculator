import sys
import argparse
from typing import Optional, Union
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from astropy import units as u
from astropy.time import Time
from rich import print as rprint
from rich_argparse import RawTextRichHelpFormatter
import plotext as pltt
from vlbiplanobs import stations
from vlbiplanobs import observation as obs
from vlbiplanobs import sources
from vlbiplanobs import freqsetups
from vlbiplanobs.gui import plots
# from vlbiplanobs import scheduler

_STATIONS = stations.Stations()
_NETWORKS = _STATIONS.get_networks_from_configfile()


def _optimal_units(value: u.Quantity, units: list[u.Unit]):
    """Given a value (with some units), returns the unit choice from all
    `units` possibilities that better suits the value.
    It is meant for the following use:
    Given 0.02*u.Jy and units = [u.kJy, u.Jy, u.mJy, u.uJy], it will
    return 20*u.mJy.
    units should have a decreasing-scale of units, and all of them
    compatible with the units in `value`.
    """
    for a_unit in units:
        if 0.8 < value.to(a_unit).value <= 800:
            return value.to(a_unit)

    # Value too high or too low
    if value.to(units[0]).value > 1:
        return value.to(units[0])

    return value.to(units[-1])


class VLBIObs(obs.Observation):
    # add __init__
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def summary(self, gui: bool = True):
        if gui:
            return self._summary_gui()

        return self._summary_tui()

    def _summary_gui(self):
        self._summary_tui()
        # TODO: for now it just calls the TUI

    def _summary_tui(self):
        """Prints the infromation for the given Observation, for testing purposes
        """
        rprint("\n[bold green]VLBI observation[/bold green]")
        rprint(f"To be conducted at {self.band.replace('cm', ' cm')} ", end='')
        if self.times is not None:
            rprint(f"from {self.times[0].strftime('%d %b %Y %H:%M')}–"
                   f"{self.times[-1].strftime('%H:%M')} UTC")
        else:
            rprint("at unspecified times.")

        if self.duration is not None:
            rprint(f"With a total duration of {_optimal_units(self.duration, [u.h, u.min, u.s]):.01f}.")

        rprint(f"\n[bold green]Setup[/bold green]")
        if None not in (self.datarate, self.bandwidth, self.subbands):
            val = _optimal_units(self.datarate, [u.Gbit/u.s, u.Mbit/u.s])
            rprint(f"\nData rate of {val.value:.0f} {val.unit.to_string('unicode')}, "
                   f"producing a total bandwidth of {_optimal_units(self.bandwidth, [u.MHz, u.GHz])}, "
                   f" divided in {self.subbands} x {int(self.bandwidth.value/self.subbands)}-"
                   f"{self.bandwidth.unit} subbands, with {self.channels} channels each, "  # type: ignore
                   f"{self.polarizations} polarization.")
        else:
            rprint("[dim]No setup (data rate, bandwidth, number of subbands) specified[/dim]")

        rprint(f"\n[bold green]Stations ({len(self.stations)})[/bold green]: "
               f"{', '.join(self.stations.station_codenames)}")
        if any([s.datarate < self.datarate for s in self.stations]):
            rprint("Note that the following stations have a reduced bandwidth:")
            for s in self.stations:
                if s.datarate < self.datarate:
                    rprint(f"    [dim]{s.codename:3}: {s.datarate.value:4.0f} "
                           f"{s.datarate.unit.to_string('unicode')} "
                           f"({int(self.subbands*s.datarate/self.datarate)} subbands)[/dim]")

        rprint("\n[bold green]Sources[/bold green]:")
        if len(self.scans) > 0:
            for ablockname, ablock in self.scans.items():
                rprint(f"    - [dim]ScanBlock[/dim] '{ablockname}'")
                rprint('      ' +
                       '\n      '.join([s.name + ' [dim](' + s.coord.to_string('hmsdms') + ')[/dim]'
                                        for s in ablock.sources()]))
        else:
            rprint("[dim]No sources defined.[/dim]")
            if self.times is not None or self.duration is not None:
                rprint("\n[bold green]Expected outcome[/bold green]:")
                val = _optimal_units(self.thermal_noise(), [u.Jy/u.beam, u.mJy/u.beam, u.uJy/u.beam])
                rprint(f"[bold]Thermal rms noise (for a +/- 45° elevation source)[/bold]: "
                       f"{val.value:.01f} {val.unit.to_string("unicode")}")
                rprint("[dim](for a +/- 45° elevation source)[/dim]")

        print('\n')

    def plot_visibility(self, gui: bool = True):
        if gui:
            return self._plot_visibility_gui()

        return self._plot_visibility_tui()

    def _plot_visibility_tui(self):
        """Show plots on the Terminal User Interface with the different sources and
        when they are visible within the observation.
        """
        if self.times is None:
            rprint("[bold green]Searching for suitable GST range "
                   "(no pre-defined observing time)[/bold green]\n")
            self.times = Time('2025-09-21', scale='utc') + np.arange(0.0, 1.005, 0.01)*u.day
            doing_gst = True
        else:
            doing_gst = False

        elevs = self.elevations()
        altaz = self.altaz()
        srcup = self.is_observable()
        srcupalways = self.is_always_observable()
        rms_noise = self.thermal_noise()
        ontarget_time = self.ontarget_time
        if doing_gst:
            # Let's put it back to the original values
            gstimes = self.gstimes
            localtimes = self.times[:]
            self.times = None

        rprint(f"[bold]The blocks are observable for:[/bold]")
        for ablockname, antbool in srcup.items():
            rprint(f"    - '{ablockname}':", end='')
            if any(srcupalways[ablockname].values()):
                ant_can = [ant for ant, b in srcupalways[ablockname].items() if b]
                if len(ant_can) >= len(self.stations):
                    rprint(f" [dim](always observable by {','.join(ant_can)})[/dim]")
                else:
                    rprint(" [dim](always observable by everyone but [/dim]", end='')
                    rprint('[dim]' + ','.join([ant for ant in self.stations.station_codenames
                                               if ant not in ant_can]) + '[/dim]')
            else:
                rprint(' [dim](nobody can observe it all the time)[/dim]')

            for anti, ant in enumerate(antbool):
                rprint(f"        {ant:3}| {''.join(['◼︎' if b else ' ' for b in antbool[ant]])}")

            rprint(f"           |-{''.join(['-' for b in antbool[ant]])}|")
            if doing_gst:
                rprint(f"           {gstimes[0].to_string(sep=':', fields=2, pad=True)} GST"
                       f"{''.join([' ' for b in antbool[ant]][:-11])}"
                       f"{(gstimes[-1] + (24*u.hourangle
                                          if np.abs(localtimes[-1].mjd - localtimes[0].mjd - 1) < 0.1
                                          else 0.0*u.hourangle)).to_string(sep=':', fields=2, pad=True)}")
            else:
                rprint(f"           {self.times.datetime[0].strftime('%H:%M'):05} UTC"
                       f"{''.join([' ' for b in antbool[ant]][:-7])}"
                       f"{self.times.datetime[-1].strftime('%H:%M'):05}")

            if doing_gst:
                when_everyone = self.when_is_observable(mandatory_stations='all', return_gst=True)[ablockname]
                if len(when_everyone) > 0:
                    rprint("\n[bold]Everyone can observe the source at: [/bold]", end='')
                    rprint(', '.join([t1.to_string(sep=':', fields=2, pad=True) + '--' +
                           t2.to_string(sep=':', fields=2, pad=True) +
                           ' GST' for t1, t2 in when_everyone]))
                else:
                    rprint("\nThe source cannot be observed by all stations at the same time.")
            else:
                when_everyone = self.when_is_observable(mandatory_stations='all')[ablockname]
                if len(when_everyone) > 0:
                    rprint("\n[bold]Everyone can observe the source at: [/bold]", end='')
                    rprint(', '.join([t1.strftime('%d %b %Y %H:%M')+'--'+t2.strftime('%H:%M') +
                           ' UTC' for t1, t2 in when_everyone]))
                else:
                    rprint("\nThe source cannot be observed by all stations at the same time.")

            min_stat = 3 if len(self.stations) > 3 else min(2, len(self.stations))
            if len(self.stations) > 2:
                rprint(f"[bold]Optimal visibility range (> {min_stat} antennas):[/bold] ", end='')
                if doing_gst:
                    rprint(', '.join([t1.to_string(sep=':', fields=2, pad=True) + '--' +
                                      (t2 + (24*u.hourangle if np.abs(t1 - t2) < 0.1*u.hourangle
                                             else 0.0*u.hourangle)).to_string(sep=':', fields=2, pad=True) +
                                      ' GST' for t1, t2
                                      in self.when_is_observable(min_stations=min_stat,
                                                                 return_gst=True)[ablockname]]))
                else:
                    rprint(', '.join([t1.strftime('%d %b %Y %H:%M')+'--'+t2.strftime('%H:%M') +
                                      ' UTC' for t1, t2
                                      in self.when_is_observable(min_stations=min_stat)[ablockname]]))

            rprint("[bold]Expected rms thermal noise for the target source: [/bold]", end='')
            for src, rms in rms_noise.items():
                if any([s.type is sources.SourceType.TARGET for s in self.sources()]):
                    if src in self.sourcenames_in_block(ablockname, sources.SourceType.TARGET):
                        val = _optimal_units(rms, [u.Jy/u.beam, u.mJy/u.beam, u.uJy/u.beam])
                        rprint(f"{src}: {val.value:.02f} {val.unit.to_string("unicode")}")
                        val = _optimal_units(ontarget_time[src], [u.h, u.min, u.s, u.ms])
                        rprint("[dim]for a total on-source time of ~ "
                               f"{val.value:.2f} {val.unit.to_string("unicode")} "
                               f"(assuming the total observing time).[/dim]")
                else:
                    if src in self.sourcenames_in_block(ablockname):
                        rprint(f"{src}: {_optimal_units(rms, [u.Jy, u.mJy, u.uJy]):.02f}")

            print('\n')

    def _plot_visibility_gui(self):
        """Show plots with the different sources and when they are visible within the
        observation
        """
        if self.scans is None:
            # rprint("No scans have been defined.")
            sys.exit(0)

        figs = plots.elevation_plot(self)
        figs.show()


def get_stations(band: str, list_networks: Optional[list[str]] = None,
                 list_stations: Optional[list[str]] = None) -> stations.Stations:
    """Returns a VLBI array including the required stations.
    Each argument is a comma-separated list of names.

    Inputs
        band : str
            The observing band. It will drop the stations that do not observe at such band.
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
    if list_networks is not None:
        try:
            networks = [_NETWORKS[n] for n in list_networks]
            for n in networks:
                for s in n.station_codenames:
                    if s not in stations and band in _STATIONS[s].bands:
                        stations.append(s)
        except KeyError:
            unknown_networks: list = [n for n in list_networks if n not in _NETWORKS]  # type: ignore
            n_networks = len(unknown_networks)  # type: ignore
            rprint(f"[bold red]The network{'s' if n_networks > 1 else ''} {', '.join(unknown_networks)}"
                   f" {'are' if n_networks > 1 else 'is'} not known.[/bold red]")
            sys.exit(1)

    if list_stations is not None:
        try:
            for s in list_stations:
                a_station = _STATIONS[s.strip()].codename
                if a_station not in stations and band in _STATIONS[a_station].bands:
                    stations.append(a_station)
        except KeyError:
            rprint(f"[bold red]The station {a_station} is not known.[/bold red]")
            sys.exit(1)

    dropped_stations = [s for s in stations if band not in _STATIONS[s].bands]
    if len(dropped_stations) > 0:
        rprint("[yellow]The following antennas were ignored "
               f"because they cannot observe at {band}: {', '.join(dropped_stations)}[/yellow]")

    return _STATIONS.filter_antennas(stations)


def main(band: str, networks: Optional[list[str]] = None,
         stations: Optional[list[str]] = None,
         src_catalog: Optional[str] = None, targets: Optional[list[str]] = None,
         start_time: Optional[Time] = None,
         duration: Optional[u.Quantity] = None, datarate: Optional[u.Quantity] = None,
         gui: bool = True):
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
        gui : bool  (default True)
            Create plots in Plotly and show them in a browser page. If False, it will only
            print some quick plots within the terminal instead.
    """
    if networks is None and stations is None:
        rprint("[bold red]You need to provide at least a VLBI network "
               "or a list of antennas that will participate in the observation.[/bold red]")
        sys.exit(1)

    if start_time is not None and start_time.scale != 'utc':
        rprint("[bold red]The start time must be in UTC[/bold red]\n"
               "[red](use 'scale' when defining the Time object)[/red]")
        sys.exit(1)

    if src_catalog is not None:
        source_catalog = sources.SourceCatalog(src_catalog)
    else:
        source_catalog = None

    src2observe: dict[str, sources.ScanBlock] = {}
    if targets is not None:
        for target in targets:
            if (source_catalog is not None) and (target in source_catalog.blocknames):
                src2observe[target] = source_catalog[target]
            else:
                a_source = sources.Source.source_from_str(target)
                src2observe[a_source.name] = sources.ScanBlock([
                    sources.Scan(a_source, duration=freqsetups.phaseref_cycle(band)
                                 if freqsetups.phaseref_cycle(band) is not None else 5*u.min)])
    elif source_catalog is None:
        # rprint("[bold red]Either a source catalog file or a list of targets must be "
        #        "provided (or both).[/bold red]")
        # TODO: no necesarily. It can provide the thermal noises for an unknown source
        # sys.exit(1)
        pass
    else:
        src2observe = source_catalog.blocks

    for target in src2observe:
        for ascan in src2observe[target]:
            if ascan.duration is None:
                match ascan.source.type:
                    case sources.SourceType.PHASECAL:
                        ascan.duration = 1.5*u.min
                    case sources.SourceType.FRINGEFINDER:
                        ascan.duration = 4.0*u.min
                    case sources.SourceType.AMPLITUDECAL:
                        ascan.duration = 4.0*u.min
                    case sources.SourceType.POLCAL:
                        ascan.duration = 5.0*u.min
                    case sources.SourceType.PULSAR:
                        ascan.duration = 5.0*u.min
                    case _:
                        ascan.duration = freqsetups.phaseref_cycle(band) - 1.5*u.min \
                            if freqsetups.phaseref_cycle(band) is not None else 5.0*u.min

    # I see a difference between different versions of Python (maybe astropy)?!!
    if duration is None:
        duration_val = None
    else:
        assert isinstance(duration, u.Quantity)
        duration_val = duration.to(u.min).value

    if datarate is None:
        if networks is not None:
            for a_network in _NETWORKS:
                if a_network in networks:
                    datarate = _NETWORKS[a_network].max_datarate(band)
                    print(f"This is a {a_network}, with {datarate})")
                    break

    o = VLBIObs(band, get_stations(band, networks, stations), scans=src2observe,
                times=start_time + np.arange(0, duration_val + 5, 10)*u.min
                if start_time is not None and duration_val is not None else None, duration=duration,
                datarate=datarate,  # subbands=, channels=, polarizations=, inttime=)
                ontarget=0.6)
    o.summary(gui)
    if targets is not None or source_catalog is not None:
        if gui:
            o.plot_visibility(gui)

        # TODO: for now I keep this so it shows the ranges in times and thermal noises
        o.plot_visibility(False)
    return o


def cli():
    usage = "%(prog)s [-h]  OPTIONS"
    description = "EVN Observation Planner"
    parser = argparse.ArgumentParser(description=description, prog="planobs", usage=usage,
                                     formatter_class=RawTextRichHelpFormatter)
    parser.add_argument('-t', '--targets', type=str, default=None, nargs='+',
                        help="Source(s) to be observed. It can be either the coordinates of the source\n"
                        "(in 'hh:mm:ss dd:mm:ss' or 'XXhXXmXXs XXdXXmXXs' format), or the name of\n"
                        "the source, if it is a known source in SIMBAD/NED/VizieR databases.\n"
                        "Or if the personal source "
                        "catalog is defined by '--input', then this\nselects the block(s) defined in the "
                        "file to use, ignoring the rest of\nsources.\nMultiple sources can be provided.")
    parser.add_argument('-i', '--input', type=str, default=None,
                        help="Input file containing the personal source catalog.\n"
                        "If provided, then '--targets' will select the block(s) "
                        "defined in\nthis file, ignoring the rest.")
    parser.add_argument('-t1', '--starttime', type=str, default=None,
                        help="Start of the observation, with the format 'YYYY-MM-DD HH:MM' "
                        "in UTC.")
    parser.add_argument('-d', '--duration', type=str, default=None,
                        help="Total duration of the observation, in hours.")
    parser.add_argument('-n', '--network', type=str, nargs='+',
                        help="The VLBI network(s) that will participate in\nthe observation. "
                        "It will take the default stations in each network.\nIf 'stations' "
                        "is provided, then it will take both the default stations\nplus "
                        "the ones given in stations. [green]See '--list-networks' to get a list.[/green]")
    parser.add_argument('-s', '--stations', type=str, nargs='+', help="List "
                        "of the antennas that will participate in the\nobservation. "
                        "You can use either antenna codenames or the standard name,\n"
                        "as given in the catalogs. [green]See '--list-antennas' to get a list.[/green]")
    parser.add_argument('-b', '--band', type=str, help="Observing band, as defined "
                        "in the catalogs as 'XXcm', with 'XX' being\nthe wavelegnth in cm. "
                        "[green]See '--list-bands' to get a list.[/green]")
    # parser.add_argument('', type=, default=, help='')
    parser.add_argument('--list-antennas', action="store_true", default=False,
                        help="Prints the list of all antennas defined in PlanObs.")
    parser.add_argument('--list-networks', action="store_true", default=False,
                        help="Prints the list of all VLBI networks defined in PlanObs.")
    parser.add_argument('--list-bands', action="store_true", default=False,
                        help="Writes the list of all observing bands defined in PlanObs.")
    # parser.add_argument('', type=, default=, help='')
    parser.add_argument('--sched', default=None, type=str,
                        help="Produces a (SCHED) .key schedule file for "
                        "the observation with the\ngiven name.")
    parser.add_argument('--fringefinders', default='2', type=str, nargs='+',
                        help="Defines the fringe finder source(s) to be scheduled "
                        "in the observation.\nIt can be either a list of source names "
                        "(as long as they\nappear in AstroGeo) "
                        "or a single number, meaning how many scans should go on\nfringe "
                        "finders, and it will automatically select the most suitable sources.")
    parser.add_argument('--polcal', action="store_true", default=False,
                        help="Requires polarization calibration for the observation.")
    parser.add_argument('--pulsar', default=None, type=str,
                        help="Sets to schedule at least a scan on a pulsar source. "
                        "If a number,\nit will select a pulsar from the personal "
                        "input source file (must be\nprovided!). If a name or list of names, "
                        "it will pick such pulsar(s).")
    parser.add_argument('--data-rate', type=float, default=None,
                        help="Maximum data rate of the observation, in Mb/s.")
    parser.add_argument('--no-gui', action="store_false", default=True,
                        help="If set, then it will not open graphical plots, but it will only\n"
                        "show the quick plots through terminal.")
    # TODO: add argument, min number of antennas possible

    args = parser.parse_args()

    if args.list_networks:
        rprint("[bold]Available VLBI networks:[/bold]")
        for network_name, network in _NETWORKS.items():
            rprint(f"[bold]{network_name}:[/bold] {network.name}")
            rprint(f"  [dim]Default antennas: {', '.join(network.station_codenames)}[/dim]")
            rprint(f"  [dim]Observes at: {', '.join(network.observing_bands)}[/dim]")

    if args.list_antennas:
        rprint("\n[bold]All available antennas:[/bold]")
        for ant in _STATIONS:
            rprint(f"     {ant.name} ({ant.codename}):  {ant.diameter} in {ant.country}")
            rprint(f"      [dim]Observes at {', '.join(ant.bands)}[/dim]")

    if args.list_bands:
        rprint("\n[bold]Available observing bands:[/bold]")
        for aband in obs.freqsetups.bands:
            rprint(f"[bold]{aband}[/bold] [dim]({obs.freqsetups.bands[aband]})[/dim]")
            rprint("[dim]  Observable with [/dim]", end='')
            rprint(f"[dim]{', '.join([nn for nn, n in _NETWORKS.items()
                                      if aband in n.observing_bands])}[/dim]")

    if args.list_antennas or args.list_bands or args.list_networks:
        sys.exit(0)

    if args.band is None:
        parser.print_help()
        rprint("\n\n[bold red]The observing band (-b/--band) and either '--network' and/or "
               "'--stations' are mandatory.[/bold red]")
        exit(1)

    if args.band not in obs.freqsetups.bands:
        rprint(f"[bold red]The provided band ({args.band}) is not available"
               "[/bold red]\n[dim]These are the available bands: "
               f"{', '.join(obs.freqsetups.bands)}.[/dim]")
        sys.exit(1)

    if args.network is None and args.stations is None:
        rprint("[bold red]You need to provide at least a VLBI network "
               "or a list of antennas that will participate in the observation.[/bold red]")
        sys.exit(1)

    if args.duration is None and args.starttime is not None:
        rprint("[bold red]If you provide a start time, you also need to provide a duration for"
               " the observation.[/bold red]")
        sys.exit(1)

    main(band=args.band, networks=args.network, stations=args.stations,
         src_catalog=args.input,
         targets=args.targets, start_time=Time(args.starttime, scale='utc')
         if args.starttime else None,
         duration=float(args.duration)*u.hour if args.duration is not None else None,
         datarate=args.data_rate*u.Mbit/u.s if args.data_rate else None, gui=args.no_gui)


if __name__ == '__main__':
    o = cli()
