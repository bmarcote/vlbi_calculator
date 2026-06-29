from __future__ import annotations
import os
import sys
import argparse
from typing import Optional, TYPE_CHECKING
from datetime import datetime as dt
from importlib.metadata import version
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from rich import print as rprint
from rich import box
from rich.table import Table
from rich.text import Text
from rich.live import Live
from rich_argparse import RawTextRichHelpFormatter

if TYPE_CHECKING:
    # Static bindings for the lazily-imported heavy names (see `_load_heavy`).
    import numpy as np
    from astropy import units as u
    from astropy.time import Time
    from astropy.coordinates import SkyCoord
    from vlbiplanobs import stations, sources, calibrators, freqsetups
    from vlbiplanobs import observation as obs
    from vlbiplanobs.cli_obs import VLBIObs, optimal_units  # noqa: F401

# Heavy dependencies (numpy, astropy, computation modules) are imported lazily via
# `_load_heavy()` so that `planobs`, `planobs -h` and `planobs -V` respond immediately.
_HEAVY_LOADED = False

# Names made available in the module globals by `_load_heavy()`.
_HEAVY_NAMES = ('np', 'u', 'Time', 'SkyCoord', 'stations', 'obs', 'sources',
                'calibrators', 'freqsetups', 'VLBIObs', 'optimal_units')


def _load_heavy() -> None:
    """Import the heavy dependencies and expose them as module globals.

    Idempotent: subsequent calls return immediately. Must be called before using
    any of the names listed in `_HEAVY_NAMES` inside this module.
    """
    global _HEAVY_LOADED
    if _HEAVY_LOADED:
        return

    import numpy as np
    from astropy import units as u
    from astropy.time import Time
    from astropy.coordinates import SkyCoord
    from vlbiplanobs import stations, sources, calibrators, freqsetups
    from vlbiplanobs import observation as obs
    from vlbiplanobs.cli_obs import VLBIObs, optimal_units

    globals().update(np=np, u=u, Time=Time, SkyCoord=SkyCoord, stations=stations,
                     obs=obs, sources=sources, calibrators=calibrators,
                     freqsetups=freqsetups, VLBIObs=VLBIObs, optimal_units=optimal_units)
    _HEAVY_LOADED = True


def __getattr__(name: str):
    """Module-level lazy attribute access (PEP 562).

    Allows `from vlbiplanobs.cli import VLBIObs` (and friends) to keep working
    while deferring the heavy imports until actually needed.
    """
    if name in _HEAVY_NAMES:
        _load_heavy()
        return globals()[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def get_stations(band: str, list_networks: Optional[list[str]] = None,
                 list_stations: Optional[list[str]] = None
                 ) -> tuple[stations.Stations, dict[str, str]]:
    """Returns a VLBI array including the required stations and any that were excluded.
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

    Returns
        tuple[Stations, dict[str, str]]
            The selected Stations object and a dict mapping excluded station codenames
            to the reason they were dropped (e.g. 'no band').
    """
    _load_heavy()
    selected = []
    no_band: dict[str, str] = {}
    if list_networks is not None:
        try:
            networks = [obs._NETWORKS[n] for n in list_networks]
            for n in networks:
                for s in n.station_codenames:
                    if s not in selected:
                        if band in obs._STATIONS[s].bands:
                            selected.append(s)
                        else:
                            no_band[s] = 'no band'
        except KeyError:
            unknown_networks: list = [n for n in list_networks if n not in obs._NETWORKS]  # type: ignore
            n_networks = len(unknown_networks)  # type: ignore
            rprint(f"[bold red]The network{'s' if n_networks > 1 else ''} {', '.join(unknown_networks)}"
                   f" {'are' if n_networks > 1 else 'is'} not known.[/bold red]")
            raise ValueError(f"Network ({unknown_networks}) not known")

    if list_stations is not None:
        try:
            for s in list_stations:
                try:
                    # Try case-sensitive lookup first
                    a_station = obs._STATIONS[s.strip()].codename
                except KeyError:
                    try:
                        # Try case-insensitive lookup by searching through all stations
                        s_upper = s.strip().upper()
                        found_station = None
                        for codename in obs._STATIONS.station_codenames:
                            if codename.upper() == s_upper:
                                found_station = obs._STATIONS[codename].codename
                                break

                        if found_station is None:
                            # Also try full station names (case insensitive)
                            for name in obs._STATIONS.station_names:
                                if name.upper() == s_upper:
                                    found_station = obs._STATIONS[name].codename
                                    break

                        if found_station is None:
                            raise KeyError(f"Station {s} not found")

                        a_station = found_station
                    except KeyError:
                        rprint(f"[bold red]The station {s} is not known.[/bold red]")
                        raise ValueError(f"Station ({s}) not known.")

                if a_station not in selected:
                    if band in obs._STATIONS[a_station].bands:
                        selected.append(a_station)
                    elif a_station not in no_band:
                        no_band[a_station] = 'no band'
        except ValueError:
            raise

    final_stations = obs._STATIONS.filter_antennas(selected)
    if not final_stations:
        rprint(f"[bold red]No antennas have been selected or none can observe at {band}.[/bold red]")
        raise ValueError("Antennas must be selected.")

    return final_stations, no_band


def _resolve_calibrators(names: list[str], target: sources.Source, band: str,
                         source_type: sources.SourceType,
                         auto_func: str = 'phasecal',
                         phase_cal_ref: Optional[sources.Source] = None) -> list[sources.Source]:
    """Resolve calibrator source names or auto-select them from the RFC catalog.

    Parameters
    ----------
    names : list[str]
        Source names to look up.  If empty, auto-selection is used.
    target : Source
        The target source (used for proximity search).
    band : str
        Observing band (e.g. '6cm').
    source_type : SourceType
        Type to assign (PHASECAL or CHECKSOURCE).
    auto_func : str
        'phasecal' or 'check' — selects which auto-selection algorithm to use.
    phase_cal_ref : Source or None
        Required when auto_func='check'; the already-selected phase calibrator.

    Returns
    -------
    list[Source]
        Resolved Source objects with the given source_type assigned.
    """
    _load_heavy()
    resolved: list[sources.Source] = []
    rfc_band = calibrators._wavelength_to_rfc_band(band)
    cat = calibrators.RFCCatalog(min_flux=0.0 * u.Jy, band=rfc_band)

    if not names:
        # Auto-select
        if auto_func == 'phasecal':
            rprint("[yellow]No phase calibrator name given — auto-selecting the best candidate. "
                   "A better source may be found manually via 'planobs phasecals'.[/yellow]")
            src = calibrators.select_phase_calibrator(target, band, catalog=cat)
        else:
            rprint("[yellow]No check source name given — auto-selecting the best candidate. "
                   "A better source may be found manually via 'planobs phasecals'.[/yellow]")
            if phase_cal_ref is None:
                phase_cal_ref = target
            src = calibrators.select_check_source(target, phase_cal_ref, band, catalog=cat)

        if src is not None:
            rprint(f"[green]  → Selected {src.name} (sep "
                   f"{target.coord.separation(src.coord).deg:.2f}°, "
                   f"unresolved {src.unresolved_flux(rfc_band):.2f} Jy)[/green]")
            resolved.append(sources.Source(
                name=src.name, coordinates=src.coord,
                source_type=source_type, other_names=[src.ivsname]))
        else:
            rprint(f"[bold red]Could not auto-select a {source_type.name.lower()} "
                   f"near {target.name}.[/bold red]")
    else:
        for sname in names:
            parsed_name, parsed_coord = sources.Source.parse_source_spec(sname)
            if parsed_name is not None and parsed_coord is not None:
                # User provided 'name/coordinates' — use explicit coordinates
                resolved.append(sources.Source(
                    parsed_name,
                    coordinates=sources.Source._parse_coord_str(parsed_coord),
                    source_type=source_type))
            else:
                lookup = parsed_name or sname
                rfc_src = cat.get_source(lookup)
                if rfc_src is not None:
                    resolved.append(sources.Source(
                        name=rfc_src.name, coordinates=rfc_src.coord,
                        source_type=source_type, other_names=[rfc_src.ivsname]))
                else:
                    try:
                        resolved.append(sources.Source.source_from_str(sname, source_type=source_type))
                    except ValueError:
                        rprint(f"[bold red]Source '{sname}' not found in RFC catalog or external "
                               f"services — skipping.[/bold red]")

    return resolved


def main(band: str, networks: Optional[list[str]] = None,
         stations: Optional[list[str]] = None, station_catalog: Optional[str] = None,
         src_catalog: Optional[str] = None, targets: Optional[list[str]] = None,
         start_time: Optional[Time] = None,
         duration: Optional[u.Quantity] = None, datarate: Optional[u.Quantity] = None,
         ontarget: float = 0.7, subbands: int = 4,
         channels: int = 64, polarizations: int = 4, inttime: Optional[u.Quantity] = None,
         phasecal_names: Optional[list[str]] = None,
         check_source_names: Optional[list[str]] = None,
         fringefinder_spec: Optional[list[str]] = None,
         polcal: bool = False) -> VLBIObs:
    """Planner for VLBI observations.

    Parameters:
    -----------

        band : str
            Observing band, as defined in the catalogs as 'XXcm', with 'XX' being
            the wavelength in cm. See .freqsetups.bands to get a list.

        stations : list of str, optional
            List of the antennas that will participate in the observation.
            You can use either antenna codenames or the standard name,
            as given in the catalogs. See .STATIONS.stations to get a list.

        station_catalog : str, optional
            Path to the file containing the list of antennas that will
            participate in the observation, if you have a local file different from the one
            distributed by PlanObs.

        src_catalog : str, optional
            Path to the toml file containing the list of sources that will be observed.
            This allows you to store a large number of sources of define them in a more detail.
            See the documentaion for help.

        targets: list[str]
            List of sources to be observed. Each entry can be:
            a) the name defining a block in the source catalog file (if provided),
            b) the coordinates of the source, in RA, DEC (J2000), as 'hh:mm:ss dd:mm:ss'
               or 'XXhXXmXXs XXdXXmXXs',
            c) the name of the source, if it is a known one so it can be found in the
               SIMBAD/NEW/VizieR databases,
            d) 'name/coordinates' to provide both a custom name and explicit coordinates
               (the coordinates override any catalog lookup).
            A mix of the previous ones can also be used for each entry.

        start_time : Time, optional
            Start of the observation, of Time class and in UTC.

        duration : astropy.units.Quantity, optional
            Total duration of the observation, as a Quantity (e.g. 1.5*u.hour).

        datarate : astropy.units.Quantity, optional
            Maximum data rate of the observation (e.g. 4*Gbit/s).

        ontarget : float (default = 0.7)
            Fraction of the total time of the observation spent on the target source.
            If multiple sources are given (e.g. already specifying scan lengths), this will be ignored.

        subbads : int (default = 4)
            Number of subbands in which the total bandwidth is split.

        channels : int (defualt = 64)
            Number of spectral channels in which each subband is divided.

        polarizations : int (defualt = 4)
            Number of polarizations recorded. It can be 1 (single pol.), 2 (dual pol.), or 4 (full stokes
            recorded: RR, LL, RL, LR).

        inttime : astropy.units.Quantity (default 2 s)
            Integration time used in the observations (e.g. time resolution on the correlated data).

    Returns:
    --------
        VLBIObs: a VLBI Observation object with all defined parameters.
    """
    _load_heavy()
    if inttime is None:
        inttime = 2.0*u.s

    if networks is None and stations is None:
        rprint("[bold red]You need to provide at least a VLBI network "
               "or a list of antennas that will participate in the observation.[/bold red]")
        raise ValueError('Either network or list of antennas must be speified.')

    # This should through an error but it will make it easier for users...
    if isinstance(networks, str):
        networks = [networks,]

    if start_time is not None and start_time.scale != 'utc':
        rprint("[bold red]The start time must be in UTC[/bold red]\n"
               "[red](use 'scale' when defining the Time object)[/red]")
        raise ValueError('Start time must be given in UTC')

    if station_catalog is not None:
        obs._STATIONS = obs.Stations(filename=station_catalog)
        obs._NETWORKS = obs.Stations.get_networks_from_configfile(stations_filename=station_catalog)

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
                try:
                    a_source = sources.Source.source_from_str(target)
                    scans_for_block: list[sources.Scan] = []

                    # Resolve phase calibrator(s)
                    if phasecal_names is not None:
                        pc_sources = _resolve_calibrators(
                            phasecal_names, a_source, band,
                            sources.SourceType.PHASECAL, auto_func='phasecal')
                        pc_dur = 1.5 * u.min
                        for pc in pc_sources:
                            scans_for_block.append(sources.Scan(pc, duration=pc_dur))

                    # Resolve check source(s)
                    if check_source_names is not None:
                        # Need a phase cal reference for geometry; use first resolved pc or target
                        pc_ref = (scans_for_block[0].source
                                  if scans_for_block else a_source)
                        cs_sources = _resolve_calibrators(
                            check_source_names, a_source, band,
                            sources.SourceType.CHECKSOURCE, auto_func='check',
                            phase_cal_ref=pc_ref)
                        cs_dur = 1.5 * u.min
                        for cs in cs_sources:
                            scans_for_block.append(
                                sources.Scan(cs, duration=cs_dur, every=4))

                    # Target scan
                    target_dur = freqsetups.phaseref_cycle(band)
                    if target_dur is not None and scans_for_block:
                        # Subtract phase-cal time from cycle to get target time
                        pc_time = sum((s.duration for s in scans_for_block
                                       if s.source.type == sources.SourceType.PHASECAL),
                                      start=0.0 * u.min)
                        target_dur = max(target_dur - pc_time, 1.0 * u.min)
                    elif target_dur is None:
                        target_dur = 5.0 * u.min
                    scans_for_block.append(
                        sources.Scan(a_source, duration=target_dur))

                    src2observe[a_source.name] = sources.ScanBlock(scans_for_block)
                except Exception:
                    raise ValueError(
                        f"Source '{target}' not found: not in the provided catalog, not in the "
                        "RFC calibrator catalog, and could not be resolved online "
                        "(SIMBAD/NED/VizieR). Check the source name or add it to your catalog.")
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
            network_station_codes = []
            for n in networks:
                if n in obs._NETWORKS:
                    for s in obs._NETWORKS[n].station_codenames:
                        if s not in network_station_codes:
                            network_station_codes.append(s)
            filtered_stations = obs._STATIONS.filter_antennas(
                network_station_codes + (stations or [])
            )
        else:
            filtered_stations = obs._STATIONS.filter_antennas(stations or [])
        if not filtered_stations:
            raise ValueError(f"No valid stations provided. Check station codes: {stations}")
        obs_networks = networks or obs.Observation.guess_network(band, filtered_stations)
        for a_network in obs._NETWORKS:
            if a_network in obs_networks:
                datarate = obs._NETWORKS[a_network].max_datarate(band)
                break
    elif isinstance(datarate, int):
        datarate = datarate*u.Mbit/u.s
        rprint("[yellow]Data rate as an int, assumed Mbit/s, but it should have had units[/yellow]")
    elif isinstance(datarate, str):
        raise ValueError("Data rate is a str! ", datarate)

    selected_stations, excluded_stations = get_stations(band, networks, stations)
    o = VLBIObs(band, selected_stations, scans=src2observe,
                times=start_time + np.arange(0, duration_val + 5, 10)*u.min
                if start_time is not None and duration_val is not None else None, duration=duration,
                datarate=datarate,
                subbands=subbands, channels=channels,
                polarizations=polarizations,
                inttime=inttime,
                ontarget=ontarget)
    o._excluded_stations = excluded_stations
    return o


def _maybe_setup_logging(logging_arg):
    """Enable file logging when the ``--logging`` flag was provided.

    Parameters
    ----------
    logging_arg : bool | str
        The parsed value of the ``--logging`` flag: False when not given,
        True when given without a path, or a string path otherwise.
    """
    if not logging_arg:
        return

    from .gui.main import setup_file_logging
    setup_file_logging(None if logging_arg is True else logging_arg)


def cli():
    """Main CLI entry point with subcommands."""
    # Handle version argument early
    if len(sys.argv) > 1 and sys.argv[1] in ('-V', '--version'):
        print(f"planobs {version('vlbiplanobs')}")
        sys.exit(0)

    # Check if this is legacy mode (no subcommand provided)
    if len(sys.argv) == 1:
        # No arguments - show subcommand help
        parser = argparse.ArgumentParser(description="EVN Observation Planner\n\n"
                       "Available modes:\n"
                       "  planobs [options]                          - Plan VLBI observations (default mode)\n"
                       "  planobs fringefinders [options]            - Find fringe finder sources\n"
                       "  planobs phasecals [options] SOURCE_NAME    - Find phase calibrator sources\n"
                       "  planobs source [options] SOURCE_NAME       - Get information about a specific source\n"
                       "  planobs server [options]                   - Start the web server\n\n"
                       "Use 'planobs <command> --help' for detailed help on each mode.",
            prog="planobs", formatter_class=RawTextRichHelpFormatter)
        parser.add_argument('-V', '--version', action='version', version=f"%(prog)s {version('vlbiplanobs')}")
        subparsers = parser.add_subparsers(dest='command', help='Available commands')
        subparsers.add_parser('observe', help='Plan VLBI observations (default mode)')
        subparsers.add_parser('fringefinders', help='Find fringe finder sources')
        subparsers.add_parser('phasecals', help='Find phase calibrator sources')
        subparsers.add_parser('server', help='Start the web server')
        parser.print_help()
        sys.exit(0)
    elif len(sys.argv) > 1 and sys.argv[1] not in ('observe', 'fringefinders', 'phasecals', 'source', 'server', 'antenna', 'ant'):
        # Legacy mode: treat as observation planning
        parser = argparse.ArgumentParser(description="EVN Observation Planner\n\n"
                       "Available modes:\n"
                       "  planobs [options]              - Plan VLBI observations (default mode)\n"
                       "  planobs fringefinders [options] - Find fringe finder sources\n"
                       "  planobs phasecals [options]     - Find phase calibrator sources\n"
                       "  planobs source [options]        - Get information about a specific source\n"
                       "  planobs server [options]        - Start the web server\n\n"
                       "Use 'planobs <command> --help' for detailed help on each mode.", prog="planobs", formatter_class=RawTextRichHelpFormatter)
        parser.add_argument('-V', '--version', action='version', version=f"%(prog)s {version('vlbiplanobs')}")
        add_observation_arguments(parser)
        add_logging_argument(parser)
        args = parser.parse_args()
        args.command = 'observe'
        _maybe_setup_logging(getattr(args, 'logging', False))
        handle_observation_command(args)
        return

    parser = argparse.ArgumentParser(
        description="EVN Observation Planner\n\n"
                   "Available modes:\n"
                   "  planobs [options]              - Plan VLBI observations (default mode)\n"
                   "  planobs fringefinders [options] - Find fringe finder sources\n"
                   "  planobs phasecals [options]     - Find phase calibrator sources\n"
                   "  planobs source [options]        - Get information about a specific source\n"
                   "  planobs server [options]        - Start the web server\n\n"
                   "Use 'planobs <command> --help' for detailed help on each mode.",
        prog="planobs", formatter_class=RawTextRichHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version=f"%(prog)s {version('vlbiplanobs')}")

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    obs_parser = subparsers.add_parser('observe', help='Plan VLBI observations (default mode)',
                                       formatter_class=RawTextRichHelpFormatter)
    add_observation_arguments(obs_parser)

    fringe_parser = subparsers.add_parser('fringefinders', help='Find fringe finder sources',
                                          formatter_class=RawTextRichHelpFormatter)
    add_fringe_finder_arguments(fringe_parser)

    phase_parser = subparsers.add_parser('phasecals', help='Find phase calibrator sources near a target',
                                         formatter_class=RawTextRichHelpFormatter)
    add_phase_cal_arguments(phase_parser)

    source_parser = subparsers.add_parser('source', help='Get information about a specific source',
                                          formatter_class=RawTextRichHelpFormatter)
    add_source_arguments(source_parser)

    server_parser = subparsers.add_parser('server', help='Start the PlanObs web server',
                                          formatter_class=RawTextRichHelpFormatter)
    add_server_arguments(server_parser)

    antenna_parser = subparsers.add_parser('antenna', aliases=['ant'],
                                           help='Get information about a specific antenna or list antennas by band',
                                           formatter_class=RawTextRichHelpFormatter)
    add_antenna_arguments(antenna_parser)

    for subparser in (obs_parser, fringe_parser, phase_parser, source_parser,
                      server_parser, antenna_parser):
        add_logging_argument(subparser)

    args = parser.parse_args()
    _maybe_setup_logging(getattr(args, 'logging', False))

    if args.command == 'observe':
        handle_observation_command(args)
    elif args.command == 'fringefinders':
        handle_fringe_finder_command(args)
    elif args.command == 'phasecals':
        handle_phase_cal_command(args)
    elif args.command == 'source':
        handle_source_command(args)
    elif args.command == 'server':
        handle_server_command(args)
    elif args.command in ('antenna', 'ant'):
        handle_antenna_command(args)


def add_observation_arguments(parser):
    """Add arguments for observation planning."""
    parser.add_argument('-t', '--targets', type=str, default=None, nargs='+',
                        help="Source(s) to be observed. Each entry can be:\n"
                        "  a) A source name (looked up in SIMBAD/NED/VizieR/RFC).\n"
                        "  b) Coordinates: 'hh:mm:ss dd:mm:ss' or 'XXhXXmXXs XXdXXmXXs'.\n"
                        "  c) 'name/coordinates' to provide both (coordinates override lookup).\n"
                        "Or if '--source-catalog' is defined, selects the block(s) in that file.\n"
                        "Multiple sources can be provided.")
    parser.add_argument('-sc', '--source-catalog', '--sc', type=str, default=None,
                        help="Input file containing the personal source catalog.\n"
                        "If provided, then '--targets' will select the block(s) "
                        "defined in\nthis file, ignoring the rest.")
    parser.add_argument('--station-catalog', type=str, default=None,
                        help="Input file containing the personal station catalog.\n"
                        "If provided, then the default catalog will not be read.")
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
    parser.add_argument('--list-antennas', action="store_true", default=False,
                        help="Prints the list of all antennas defined in PlanObs.")
    parser.add_argument('--list-networks', action="store_true", default=False,
                        help="Prints the list of all VLBI networks defined in PlanObs.")
    parser.add_argument('--list-bands', action="store_true", default=False,
                        help="Writes the list of all observing bands defined in PlanObs.")
    parser.add_argument('--sched', default=None, type=str,
                        help="Produces a (SCHED) .key schedule file for "
                        "the observation with the\ngiven name.")
    parser.add_argument('--fringefinders', default='2', type=str, nargs='+',
                        help="Defines the fringe finder source(s) to be scheduled "
                        "in the observation.\nIt can be either a list of source names "
                        "(as long as they\nappear in AstroGeo), "
                        "'name/coordinates' to provide both,\n"
                        "or a single number, meaning how many scans should go on\nfringe "
                        "finders, and it will automatically select the most suitable sources.")
    parser.add_argument('--polcal', action="store_true", default=False,
                        help="Requires polarization calibration for the observation.")
    parser.add_argument('--phasecal', default=None, type=str, nargs='*',
                        help="Phase calibrator source(s) for the target. If no names are given\n"
                        "(just --phasecal), the best candidate is picked automatically.\n"
                        "One or more source names or 'name/coordinates' can be provided.")
    parser.add_argument('--check-source', default=None, type=str, nargs='*',
                        help="Check source(s) for the target. If no names are given\n"
                        "(just --check-source), the best candidate is picked automatically.\n"
                        "One or more source names or 'name/coordinates' can be provided.")
    parser.add_argument('--pulsar', default=None, type=str,
                        help="Sets to schedule at least a scan on a pulsar source. "
                        "If a number,\nit will select a pulsar from the personal "
                        "input source file (must be\nprovided!). If a name, 'name/coordinates', "
                        "or coordinates, it will\nresolve accordingly.")
    parser.add_argument('--data-rate', type=float, default=None,
                        help="Maximum data rate of the observation, in Mb/s.")
    parser.add_argument('--gui', action="store_true", default=False,
                        help="If set, then it will not open graphical plots, but it will only\n"
                        "show the quick plots through terminal.")
    parser.add_argument('--no-tui', action="store_false", default=True,
                        help="If set, then it will not show all the output in the terminal as default.")
    parser.add_argument('--debug', action="store_true", default=False,
                        help="If set, shows some debuging messages.")


def add_fringe_finder_arguments(parser):
    """Add arguments for fringe finder search."""
    # Import the main_fringe function to get its argument parser
    # We'll recreate the arguments here
    parser.add_argument('-n', '--network', type=str, nargs='+',
                        help="The VLBI network(s) that will participate in\nthe observation. "
                        "It will take the default stations in each network.\nIf 'stations' "
                        "is provided, then it will take both the default stations\nplus "
                        "the ones given in stations. [green]See '--list-networks' to get a list.[/green]")
    parser.add_argument('-s', '--stations', type=str, nargs='+',
                        help="List of antenna codenames or names that will participate in the observation.")
    parser.add_argument('-t', '--starttime', type=str, required=True,
                        help="Start of the observation in format 'YYYY-MM-DD HH:MM' (UTC).")
    parser.add_argument('-d', '--duration', type=float, required=True, help="Duration of the observation in hours.")
    parser.add_argument('--min-flux', type=float, default=0.5,
                        help="Minimum unresolved flux threshold in Jy (default: 0.5).")
    parser.add_argument('--min-elevation', type=float, default=20.0,
                        help="Minimum elevation in degrees (default: 20).")
    parser.add_argument('-l', '--max-lines', type=int, default=20,
                        help="Maximum number of sources to return (default: 20).")
    parser.add_argument('--require-all', action='store_true', default=False,
                        help="Require source to be visible by ALL stations (default: False).")
    parser.add_argument('-b', '--band', type=str, default=None,
                        help="Observing band for flux display (e.g., '18cm', '6cm'). If not provided, shows flux for all available bands.")
    parser.add_argument('--station-catalog', type=str, default=None, help="Path to custom station catalog file.")
    parser.add_argument('--json', action='store_true', default=False,
                        help="Output results in JSON format instead of a table.")


def add_phase_cal_arguments(parser):
    """Add arguments for phase calibrator search."""
    parser.add_argument('-t', '--target', type=str, required=True,
                        help="Target source name (J2000 or IVS name from RFC catalog;\n"
                        "or a block/source name from '--source-catalog' if provided).")
    parser.add_argument('-sc', '--source-catalog', '--sc', type=str, default=None,
                        help="Input file containing the personal source catalog.\n"
                        "If provided, then '--target' will first be looked up in\nthis file "
                        "(block or source name), as in the observe mode.")
    parser.add_argument('--max-separation', type=float, default=5.0,
                        help="Maximum angular separation in degrees (default: 5.0).")
    parser.add_argument('--min-flux', type=float, default=0.1,
                        help="Minimum unresolved flux threshold in Jy (default: 0.1).")
    parser.add_argument('-n', '--n-sources', type=int, default=None,
                        help="Maximum number of sources to return (default: all).")
    parser.add_argument('-b', '--band', type=str, default=None,
                        help="Observing band for flux display (e.g., '18cm', '6cm'). If not provided, shows flux for all available bands.")
    parser.add_argument('--catalog-file', type=str, default=None,
                        help="Path to custom RFC catalog file.")
    parser.add_argument('--json', action='store_true', default=False,
                        help="Output results in JSON format instead of a table.")


def add_source_arguments(parser):
    """Add arguments for source information lookup."""
    parser.add_argument('source_name', type=str,
                        help="Name of the source to get information about.")
    parser.add_argument('--no-networks', action='store_true', default=False,
                        help="Skip the network observability table calculation and display.")
    parser.add_argument('--gst', action='store_true', default=False,
                        help="Show GST time ranges (HH:MM-HH:MM) when the source is visible by >3 antennas per network.")


def add_server_arguments(parser):
    """Add arguments for server."""
    parser.add_argument('--host', type=str, default='127.0.0.1', help="Host address (default: 127.0.0.1)")
    parser.add_argument('--port', type=int, default=8050, help="Port number (default: 8050)")
    parser.add_argument('--debug', action='store_true', default=False, help="Enable debug mode")


def add_logging_argument(parser):
    """Add the optional ``--logging`` flag that enables writing a log file.

    No log file is created by default. When the flag is given without a value,
    a default path is used ('/var/log/planobs.log' if writable, otherwise
    '~/log-planobs.log'). An explicit path may also be provided.
    """
    parser.add_argument('--logging', nargs='?', const=True, default=False, metavar='LOGFILE',
                        help="Enable logging to a file. Optionally provide a path;\n"
                        "otherwise '/var/log/planobs.log' (if writable) or\n"
                        "'~/log-planobs.log' is used. Disabled by default.")


def handle_observation_command(args):
    """Handle the observation planning command."""
    _load_heavy()
    t0 = dt.now() if args.debug else None

    if args.list_networks:
        rprint("[bold]Available VLBI networks:[/bold]")
        for network_name, network in obs._NETWORKS.items():
            rprint(f"[bold]{network_name}:[/bold] {network.name}")
            rprint(f"  [dim]Default antennas: {', '.join(network.station_codenames)}[/dim]")
            rprint(f"  [dim]Observes at: {', '.join(network.observing_bands)}[/dim]")

    if args.list_antennas:
        rprint("\n[bold]All available antennas:[/bold]")
        for ant in obs._STATIONS:
            rprint(f"     {ant.name} ({ant.codename}):  {ant.diameter} in {ant.country}")
            rprint(f"      [dim]Observes at {', '.join(ant.bands)}[/dim]")

    if args.list_bands:
        rprint("\n[bold]Available observing bands:[/bold]")
        for aband in obs.freqsetups.bands:
            rprint(f"[bold]{aband}[/bold] [dim]({obs.freqsetups.bands[aband]})[/dim]")
            rprint("[dim]  Observable with [/dim]", end='')
            rprint("[dim]" +
                   ', '.join([nn for nn, n in obs._NETWORKS.items() if aband in n.observing_bands]) +
                   "[/dim]")

    if args.list_antennas or args.list_bands or args.list_networks:
        sys.exit(0)

    if args.band is None:
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

    if (not args.gui) and (not args.no_tui):
        rprint("[bold yellow]Note that you supressed both GUI and TUI. "
               "No output will be provided.[/bold yellow]")

    # Resolve phasecal / check-source arguments (None = not requested, [] = auto-select)
    phasecal_arg = getattr(args, 'phasecal', None)
    check_source_arg = getattr(args, 'check_source', None)
    fringefinder_arg = getattr(args, 'fringefinders', ['2'])
    polcal_arg = getattr(args, 'polcal', False)

    try:
        o = main(band=args.band, networks=args.network, stations=args.stations,
            src_catalog=args.source_catalog, station_catalog=args.station_catalog,
            targets=args.targets, start_time=Time(args.starttime, scale='utc')
            if args.starttime else None,
            duration=float(args.duration)*u.hour if args.duration is not None else None,
            datarate=args.data_rate*u.Mbit/u.s if args.data_rate else None,
            phasecal_names=phasecal_arg,
            check_source_names=check_source_arg,
            fringefinder_spec=fringefinder_arg,
            polcal=polcal_arg)
    except ValueError as e:
        rprint(f"[bold red]Error: {e}[/bold red]")
        sys.exit(1)

    o.summary(args.gui, args.no_tui)
    if args.targets is not None or args.source_catalog is not None:
        o.plot_visibility(args.gui, args.no_tui)

    if args.sched is not None:
        key_filename = args.sched if args.sched.endswith('.key') else f"{args.sched}.key"
        from vlbiplanobs.scheduler import ObservationScheduler
        scheduler = ObservationScheduler(
            o, fringefinder_spec=fringefinder_arg, polcal=polcal_arg)
        scheduler.schedule()
        key_content = scheduler.generate_key_file(
            experiment_code=args.sched.replace('.key', '').upper())
        with open(key_filename, 'w') as f:
            f.write(key_content)
        rprint(f"[green]Schedule file written to: {key_filename}[/green]")

    if args.debug:
        print(f"Execution time: {(dt.now() - t0).total_seconds()} s")


def handle_fringe_finder_command(args):
    """Handle the fringe finder command."""
    _load_heavy()
    if args.network is None and args.stations is None:
        rprint("[bold red]You need to provide at least a VLBI network "
               "or a list of antennas that will participate in the observation.[/bold red]")
        sys.exit(1)

    original_argv = sys.argv.copy()
    sys.argv = ['planobs_fringefinder']
    if args.network is not None:
        sys.argv.extend(['-n'] + args.network)
    if args.stations is not None:
        sys.argv.extend(['-s'] + args.stations)
    sys.argv.extend(['-t', args.starttime, '-d', str(args.duration)])
    if args.min_flux != 0.5:
        sys.argv.extend(['--min-flux', str(args.min_flux)])
    if args.min_elevation != 20.0:
        sys.argv.extend(['--min-elevation', str(args.min_elevation)])
    if args.max_lines != 20:
        sys.argv.extend(['-l', str(args.max_lines)])
    if args.require_all:
        sys.argv.append('--require-all')
    if args.band is not None:
        sys.argv.extend(['-b', args.band])
    if args.station_catalog is not None:
        sys.argv.extend(['--station-catalog', args.station_catalog])
    if args.json:
        sys.argv.append('--json')
    try:
        calibrators.main_fringe()
    finally:
        sys.argv = original_argv


def handle_phase_cal_command(args):
    """Handle the phase calibrator command."""
    _load_heavy()
    original_argv = sys.argv.copy()
    sys.argv = ['planobs_phasecal', '-t', args.target]
    if args.source_catalog is not None:
        sys.argv.extend(['-sc', args.source_catalog])
    if args.max_separation != 5.0:
        sys.argv.extend(['--max-separation', str(args.max_separation)])
    if args.min_flux != 0.1:
        sys.argv.extend(['--min-flux', str(args.min_flux)])
    if args.n_sources is not None:
        sys.argv.extend(['-n', str(args.n_sources)])
    if args.band is not None:
        sys.argv.extend(['-b', args.band])
    if args.catalog_file is not None:
        sys.argv.extend(['--catalog-file', args.catalog_file])
    if args.json:
        sys.argv.append('--json')
    try:
        calibrators.main_phasecal()
    finally:
        sys.argv = original_argv


def _format_gst_ranges(gst_pairs: list[tuple]) -> str:
    """Format a list of GST (Longitude) start/end pairs as 'HH:MM-HH:MM, ...' string."""
    parts = []
    for pair in gst_pairs:
        t0, t1 = pair
        h0, m0 = int(t0.hour), int((t0.hour % 1) * 60)
        h1, m1 = int(t1.hour), int((t1.hour % 1) * 60)
        parts.append(f"{h0:02d}:{m0:02d}-{h1:02d}:{m1:02d}")
    return ", ".join(parts)


def _check_obs_worker(args):
    """Module-level worker for ProcessPoolExecutor: checks if a network can observe a source.

    Must be defined at module level (not as a closure) to be picklable.
    Accepts (net_key, band, source_name, ra_deg, dec_deg[, return_gst]).
    Returns (net_key, band, is_observable, gst_ranges_str | None).
    """
    _load_heavy()
    return_gst = args[5] if len(args) > 5 else False
    net_key, band, source_name, ra_deg, dec_deg = args[:5]
    try:
        coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
        src = sources.Source(name=source_name, coordinates=coord)
        network = obs._NETWORKS[net_key]
        observation = obs.Observation(
            band=band, stations=network, times=None, duration=1 * u.hour,
            datarate=network.max_datarate(band),
            scans={src.name: sources.ScanBlock([sources.Scan(src, duration=5 * u.min)])}
        )
        observable_times = observation.when_is_observable(min_stations=3)
        is_obs = len(observable_times[src.name]) > 0
        if return_gst and is_obs:
            gst_times = observation.when_is_observable(min_stations=3, return_gst=True)
            gst_str = _format_gst_ranges(gst_times.get(src.name, []))
            return net_key, band, True, gst_str
        return net_key, band, is_obs, "" if return_gst else None
    except Exception:
        return net_key, band, False, "" if return_gst else None


def _build_obs_table(sorted_bands: list[str],
                     observability: dict[tuple[str, str], bool],
                     pending: set[tuple[str, str]],
                     show_gst: bool = False,
                     gst_ranges: dict[tuple[str, str], str] | None = None) -> Table:
    """Build the network observability Rich Table from current results.

    Cells that are still being computed are shown as '…'.
    When show_gst is True, appends a 'GST' column per network showing GST ranges where >3 antennas observe.
    """
    table = Table(box=None, show_header=True, padding=0, pad_edge=False)
    table.add_column("Network", style="cyan", justify="left", width=10, no_wrap=True)
    for band in sorted_bands:
        table.add_column(band.replace('cm', ''), justify="center", width=6, no_wrap=True)
    if show_gst:
        table.add_column("GST", style="white", justify="left", no_wrap=False)
    for net_key, network in obs._NETWORKS.items():
        row: list[Text | str] = [net_key]
        for band in sorted_bands:
            if band in network.observing_bands:
                if (net_key, band) in pending:
                    row.append(Text("  … ", style="dim", justify="center"))
                elif observability.get((net_key, band), False):
                    row.append(Text(" ✓ ", style="bold on #90EE90", justify="center"))
                else:
                    row.append("   ")
            else:
                row.append("   ")
        if show_gst:
            net_gst_parts: list[str] = []
            any_pending = False
            for band in sorted_bands:
                if band in network.observing_bands:
                    if (net_key, band) in pending:
                        any_pending = True
                    elif gst_ranges and (net_key, band) in gst_ranges and gst_ranges[(net_key, band)]:
                        net_gst_parts.append(gst_ranges[(net_key, band)])
            if any_pending:
                row.append(Text(" … ", style="dim"))
            elif net_gst_parts:
                row.append(", ".join(dict.fromkeys(net_gst_parts)))
            else:
                row.append("")
        table.add_row(*row)
    return table


def _show_observability_table(source: sources.Source, show_gst: bool = False) -> None:
    """Compute and display the network observability table with live progressive updates.

    Tries ProcessPoolExecutor first (true parallelism, bypasses GIL) and falls back
    to ThreadPoolExecutor if process spawning fails. Results are displayed as they
    arrive via Rich Live, so the table fills in progressively rather than all at once.

    Inputs
        - source : sources.Source — the source to check observability for.
        - show_gst : bool — if True, appends a GST column with time ranges where >3 antennas observe.
    """
    all_bands: set[str] = set()
    for network in obs._NETWORKS.values():
        all_bands.update(network.observing_bands)

    def wavelength_key(band: str) -> float:
        return -float(band[:-2]) if band.endswith('cm') else 0.0

    sorted_bands = sorted(all_bands, key=wavelength_key)

    tasks = [
        (net_key, band, source.name, source.coord.ra.deg, source.coord.dec.deg, show_gst)
        for net_key, network in obs._NETWORKS.items()
        for band in sorted_bands
        if band in network.observing_bands
    ]

    observability: dict[tuple[str, str], bool] = {}
    gst_ranges: dict[tuple[str, str], str] = {}
    pending: set[tuple[str, str]] = {(t[0], t[1]) for t in tasks}
    n_workers = min(len(tasks), (os.cpu_count() or 4))

    def _run(executor_cls, worker_tasks: list) -> None:
        with executor_cls(max_workers=n_workers) as executor:
            futures = {executor.submit(_check_obs_worker, task): (task[0], task[1])
                       for task in worker_tasks}
            for future in as_completed(futures):
                try:
                    net_key, band, result, gst_str = future.result()
                except Exception:
                    net_key, band = futures[future]
                    result = False
                    gst_str = None
                observability[(net_key, band)] = result
                if gst_str is not None:
                    gst_ranges[(net_key, band)] = gst_str
                pending.discard((net_key, band))
                live.update(_build_obs_table(sorted_bands, observability, pending,
                                             show_gst=show_gst, gst_ranges=gst_ranges))

    with Live(_build_obs_table(sorted_bands, observability, pending, show_gst=show_gst, gst_ranges=gst_ranges),
              refresh_per_second=8) as live:
        try:
            _run(ProcessPoolExecutor, tasks)
        except Exception:
            # ProcessPoolExecutor failed (e.g. pickling issue on this platform);
            # retry any remaining tasks with threads.
            remaining = [t for t in tasks if (t[0], t[1]) not in observability]
            if remaining:
                _run(ThreadPoolExecutor, remaining)


def handle_source_command(args):
    """Handle the source information command."""
    _load_heavy()
    try:
        source_obj = None
        calibrator = calibrators.RFCCatalog(min_flux=0.0).get_source(args.source_name)
        if calibrator:
            rprint(f"[bold]{calibrator.name}[/bold] (also known as {calibrator.ivsname})")
            rprint(f"[bold]Coordinates:[/bold] {calibrator.coord.to_string('hmsdms')}")
            rprint("\n[bold green]AstroGeo Information[/bold green]")
            rprint(f"[bold]Number of Observations:[/bold] {calibrator.n_observations}")
            band_names = {'s': '18/21cm', 'c': '13/6/5cm', 'x': '3.6cm', 'u': '2cm', 'k': '1.3/0.7cm'}
            table = Table(box=box.SIMPLE)
            table.add_column("Band", style="cyan", justify="center")
            table.add_column("Wavelength", style="white", justify="center")
            table.add_column("Total Flux (Jy)", style="green", justify="right")
            table.add_column("Unresolved Flux (Jy)", style="yellow", justify="right")
            for band in ('s', 'c', 'x', 'u', 'k'):
                total_flux, unresolved_flux = calibrator.get_flux_at_band(band)
                if total_flux > 0 or unresolved_flux > 0:
                    table.add_row(band.upper(), band_names[band], f"{total_flux:.2f}",
                                  f"{unresolved_flux:.2f}")
            rprint(table)
            rprint(f"\n[bold][link={calibrator.get_astrogeo_link()}]AstroGeo Link[/bold]")
            # Create source object from calibrator
            source_obj = sources.Source(
                name=calibrator.name,
                coordinates=calibrator.coord,
                other_names=[calibrator.ivsname]
            )
        else:
            try:
                source_obj = sources.Source.source_from_str(args.source_name)
                rprint("\n[bold green]Source Information[/bold green]")
                rprint(f"[bold]Name:[/bold] {source_obj.name}")
                rprint(f"[bold]Coordinates:[/bold] {source_obj.coord.to_string('hmsdms')}")
            except Exception as e:
                rprint(f"[bold red]Source '{args.source_name}' not recognized[/bold red]")
                rprint(f"[red]Error: Not found in the RFC Catalog and {e} [/red]")
                sys.exit(1)

        # Show observability table for all sources
        if source_obj and not args.no_networks:
            rprint("\n[bold green]Observable by (bands in cm)[/bold green]")
            _show_observability_table(source_obj, show_gst=args.gst)
    except Exception as e:
        rprint(f"[bold red]Error:[/bold red] {e}")
        sys.exit(1)


def handle_server_command(args):
    """Handle the server command."""
    # Lazy import: pulls in dash and the whole GUI stack.
    from vlbiplanobs.gui.main import main as gui_main
    gui_main(debug=args.debug, host=args.host, port=args.port)


def _render_horizontal_band_table(bands_to_show: list[str], ant) -> None:
    """Renders antenna band/SEFD info as a horizontal table, wrapping to terminal width.

    First column shows row labels ('Band (cm)', 'SEFD (Jy)'); each subsequent column is one band.
    Splits into multiple sub-tables if total width exceeds terminal width.
    - bands_to_show: list of band strings to display (e.g. ['6cm', '18cm'])
    - ant: Station object with sefd(band) method
    """
    from rich.console import Console as _Console
    term_width = _Console().width

    col_data: list[tuple[str, str]] = []
    for b in bands_to_show:
        sefd_val = ant.sefd(b)
        sefd_str = f"{sefd_val.value:.0f}" if hasattr(sefd_val, 'value') else str(sefd_val)
        col_data.append((b, sefd_str))

    # Column width = widest of (band name, sefd value); +2 for Rich's cell padding (1 each side)
    # +1 for the column separator, giving the true rendered width per data column
    LABEL_CONTENT = max(len("Band (cm)"), len("SEFD (Jy)"))
    CELL_OVERHEAD = 3  # 1 left-pad + 1 right-pad + 1 separator (box.SIMPLE)
    label_width = LABEL_CONTENT + CELL_OVERHEAD
    col_widths = [max(len(b), len(s)) + CELL_OVERHEAD for b, s in col_data]

    # Split columns into chunks that each fit within terminal width
    chunks: list[list[tuple[str, str]]] = []
    current: list[tuple[str, str]] = []
    used = label_width
    for (b, s), w in zip(col_data, col_widths):
        if current and used + w > term_width:
            chunks.append(current)
            current = [(b, s)]
            used = label_width + w
        else:
            current.append((b, s))
            used += w
    if current:
        chunks.append(current)

    for chunk in chunks:
        table = Table(box=box.SIMPLE, show_header=False, padding=(0, 1))
        table.add_column("", style="bold magenta", no_wrap=True, min_width=LABEL_CONTENT)
        for b, s in chunk:
            cw = max(len(b), len(s))
            table.add_column(b, justify="right", no_wrap=True, min_width=cw)
        table.add_row("Band (cm)", *[f"[bold cyan]{b.replace('cm', '')}[/bold cyan]" for b, _ in chunk])
        table.add_row("SEFD (Jy)", *[s for _, s in chunk])
        rprint(table)


def add_antenna_arguments(parser):
    """Add arguments for the antenna info/listing command."""
    parser.add_argument('antenna_name', type=str, nargs='?', default=None,
                        help="Name, short name, or codename of the antenna to look up. "
                        "If omitted, lists all antennas (or all antennas at the given --band).")
    parser.add_argument('-b', '--band', type=str, default=None,
                        help="Filter by observing band (e.g. '18cm', '6cm'). "
                        "If no antenna is given, lists all antennas that observe at this band.")


def _normalize_band(band: str) -> str:
    """Normalize band string to include 'cm' suffix if missing.
    Returns the normalized band string (e.g. '18' -> '18cm', '18cm' -> '18cm').
    Only appends 'cm' if the string contains no unit suffix (no letters).
    """
    band = band.strip()
    if band.replace('.', '').isdigit():
        band = f"{band}cm"
    return band


def _find_antenna(name: str) -> Optional[object]:
    """Search obs._STATIONS for an antenna matching name, fullname, or codename (case-insensitive).
    Returns the matching Station object or None if not found.
    """
    name_lower = name.lower()
    for ant in obs._STATIONS:
        if (ant.name.lower() == name_lower or ant.fullname.lower() == name_lower
                or ant.codename.lower() == name_lower):
            return ant
    return None


def handle_antenna_command(args):
    """Handle the 'antenna'/'ant' subcommand.

    Behavior:
    - antenna_name + band: show info for that antenna filtered to the given band.
    - antenna_name only: show full info for that antenna.
    - band only: list all antennas that can observe at that band.
    - neither: list all antennas (same as --list-antennas).
    """
    _load_heavy()
    band = _normalize_band(args.band) if args.band else None
    if band is not None and band not in obs.freqsetups.bands:
        rprint(f"[bold red]Band '{band}' is not recognized.[/bold red] "
               f"Available bands: {', '.join(obs.freqsetups.bands)}")
        sys.exit(1)

    if args.antenna_name is None:
        if band is None:
            # Same as --list-antennas
            rprint("\n[bold]All available antennas:[/bold]")
            for ant in obs._STATIONS:
                rprint(f"     [bold]{ant.name}[/bold] ([cyan]{ant.codename}[/cyan]):  "
                       f"{ant.diameter} in {ant.country}")
                rprint(f"      [dim]Observes at {', '.join(ant.bands)}[/dim]")
        else:
            matching = [ant for ant in obs._STATIONS if ant.has_band(band)]
            if not matching:
                rprint(f"[bold red]No antennas found that observe at {band}.[/bold red]")
                sys.exit(1)
            rprint(f"\n[bold]Antennas that observe at [cyan]{band}[/cyan]:[/bold]")
            table = Table(box=box.SIMPLE, show_header=True, header_style="bold magenta")
            table.add_column("Antenna", style="bold")
            table.add_column("Codename", style="cyan")
            table.add_column("Diameter")
            table.add_column("Country")
            table.add_column(f"SEFD at {band} (Jy)", justify="right")
            for ant in sorted(matching, key=lambda a: a.name):
                sefd_val = ant.sefd(band)
                sefd_str = f"{sefd_val.value:.0f}" if hasattr(sefd_val, 'value') else str(sefd_val)
                table.add_row(ant.name, ant.codename, ant.diameter, ant.country, sefd_str)
            rprint(table)
        return

    ant = _find_antenna(args.antenna_name)
    if ant is None:
        rprint(f"[bold red]Antenna '{args.antenna_name}' not found.[/bold red] "
               "Run [bold]planobs --list-antennas[/bold] to see the available antennas.")
        sys.exit(1)

    rprint(f"\n[bold underline]{ant.fullname}[/bold underline]")
    rprint(f"  [bold]Short name:[/bold]  {ant.name}")
    rprint(f"  [bold]Codename:[/bold]    [cyan]{ant.codename}[/cyan]")
    rprint(f"  [bold]Diameter:[/bold]    {ant.diameter}")
    rprint(f"  [bold]Country:[/bold]     {ant.country}")

    if band is not None:
        if not ant.has_band(band):
            rprint(f"\n[bold red]{ant.name} does not observe at {band}.[/bold red] "
                   f"It observes at: {', '.join(ant.bands)}")
            sys.exit(1)
        bands_to_show = [band]
    else:
        bands_to_show = sorted(ant.bands, key=lambda b: float(b.replace('cm', '')))

    # Print SEFD table (horizontal layout, wraps at terminal width)
    rprint("")
    _render_horizontal_band_table(bands_to_show, ant)


if __name__ == '__main__':
    o = cli()
