from __future__ import annotations
import os
import re
import sys
import argparse
from pathlib import Path
from typing import Optional, TYPE_CHECKING
from datetime import datetime as dt
from importlib.metadata import version
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from rich import print as rprint
from rich import box
from rich.console import Console
from rich.table import Table
from rich.text import Text
from rich.live import Live
from rich.markup import escape
from rich_argparse import RawTextRichHelpFormatter
from vlbiplanobs.cli_options import add_deprecated_alias, add_removed_option

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
_HEAVY_NAMES = ('np', 'u', 'Time', 'SkyCoord', 'stations', 'obs', 'sources',
                'calibrators', 'freqsetups', 'VLBIObs', 'optimal_units')

# Server-safety bounds for the observation duration accepted by `main` (also used by the web GUI).
MAX_DURATION_HOURS = 96.0

# Valid SCHED experiment codes (derived from the --sched filename).
_EXPERIMENT_CODE_PATTERN = r'[A-Za-z0-9_]+'

# Help text shared by every top-level `planobs` parser (rich markup: literal brackets are escaped).
_MODES_HELP = ("EVN Observation Planner\n\n"
               "Available modes:\n"
               "  planobs \\[options]                    - Plan VLBI observations (default mode)\n"
               "  planobs fringefinders \\[options]      - Find fringe finder sources\n"
               "  planobs phasecals \\[options] TARGET   - Find phase calibrator sources\n"
               "  planobs source \\[options] SOURCE      - Get information about a specific source\n"
               "  planobs antenna \\[options] \\[ANTENNA]  - Antenna information (alias: ant)\n"
               "  planobs server \\[options]             - Start the web server\n\n"
               "Global options:\n"
               "  --get-key-template FILENAME          - Copy the bundled SCHED key template\n\n"
               "Use 'planobs <command> --help' for detailed help on each mode.")


class PagerHelpAction(argparse.Action):
    """Show help in a terminal pager if output is a TTY.

    Falls back to printing directly when the terminal has no pager or is not
    interactive (e.g. piped to a file).
    """

    def __call__(self, parser, namespace, values, option_string=None):
        help_text = Text.from_ansi(parser.format_help())
        console = Console()
        if console.is_terminal:
            lines = help_text.split('\n')
            page_height = max(console.height - 2, 1)
            for page_start in range(0, len(lines), page_height):
                page = Text('\n').join(lines[page_start:page_start + page_height])
                console.print(page)
                if page_start + page_height < len(lines):
                    try:
                        response = console.input('[dim]-- More -- (Enter to continue, q to quit) [/dim]')
                    except (EOFError, KeyboardInterrupt):
                        break
                    if response.strip().lower() == 'q':
                        break
        else:
            console.print(help_text, end='')
        parser.exit()


def _add_pager_help(parser: argparse.ArgumentParser) -> None:
    """Add a paged ``-h/--help`` action to a parser."""
    parser.add_argument('-h', '--help', action=PagerHelpAction, nargs=0,
                        help='Show this help message and exit')


def _parser_with_pager(*args, **kwargs) -> argparse.ArgumentParser:
    """Create an ArgumentParser whose ``-h/--help`` is shown in a pager."""
    parser = argparse.ArgumentParser(*args, add_help=False, **kwargs)
    _add_pager_help(parser)
    return parser


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


def copy_key_template(destination: str) -> None:
    """Copy the bundled SCHED .key template to a user-specified path.

    Parameters
    ----------
    destination : str
        Path where the template should be written.
    """
    from importlib import resources
    from pathlib import Path
    template_path = resources.files('vlbiplanobs.data').joinpath('key_file.key.template')
    Path(destination).write_text(Path(template_path).read_text(), encoding='utf-8')
    rprint(f"[green]Template written to: {escape(destination)}[/green]")


def _maybe_handle_get_key_template() -> bool:
    """Handle ``--get-key-template`` before full argument parsing.

    Returns True if the option was present and the template was copied, in which
    case the program should exit.
    """
    pre_parser = argparse.ArgumentParser(add_help=False)
    pre_parser.add_argument('--get-key-template', type=str, dest='get_key_template')
    pre_args, _ = pre_parser.parse_known_args()
    if pre_args.get_key_template:
        copy_key_template(pre_args.get_key_template)
        return True
    return False


def __getattr__(name: str):
    """Module-level lazy attribute access (PEP 562).

    Allows `from vlbiplanobs.cli import VLBIObs` (and friends) to keep working
    while deferring the heavy imports until actually needed.

    Parameters
    ----------
    name : str
        Attribute name to access.

    Returns
    -------
    object
        The requested attribute.

    Raises
    ------
    AttributeError
        If the attribute is not found.
    """
    if name in _HEAVY_NAMES:
        _load_heavy()
        return globals()[name]

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def get_stations(band: str, list_networks: Optional[list[str]] = None,
                 list_stations: Optional[list[str]] = None) -> tuple[stations.Stations, dict[str, str]]:
    """Returns a VLBI array including the required stations and any that were excluded.

    Parameters
    ----------
    band : str
        The observing band. It will drop the stations that do not observe at such band.
    list_networks : list[str] or None
        If you want to pick the default antennas participating in one of the
        known VLBI networks, then you can include directly the network name here.
    list_stations : list[str] or None
        If you want a particular list of stations, or adding some that are not
        in the default network, then you can quote them here, using either the
        station code names or their names.

    Returns
    -------
    tuple[Stations, dict[str, str]]
        The selected Stations object and a dict mapping excluded station codenames
        to the reason they were dropped (e.g. 'no band').
    """
    _load_heavy()
    selected: list[str] = []
    no_band: dict[str, str] = {}
    if list_networks:
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
            unknown_networks: list = [n for n in list_networks if n not in obs._NETWORKS]
            n_networks: int = len(unknown_networks)
            rprint(f"[bold red]The network{'s' if n_networks > 1 else ''} {escape(', '.join(unknown_networks))}"
                   f" {'are' if n_networks > 1 else 'is'} not known.[/bold red]")
            raise ValueError(f"Network ({unknown_networks}) not known")

    if list_stations:
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
                    rprint(f"[bold red]The station {escape(s)} is not known.[/bold red]")
                    raise ValueError(f"Station ({s}) not known.")

            if a_station not in selected:
                if band in obs._STATIONS[a_station].bands:
                    selected.append(a_station)
                elif a_station not in no_band:
                    no_band[a_station] = 'no band'

    final_stations = obs._STATIONS.filter_antennas(selected)
    if not final_stations:
        rprint(f"[bold red]No antennas have been selected or none can observe at {band}.[/bold red]")
        raise ValueError("Antennas must be selected.")

    return final_stations, no_band


def _rfc_catalog_for_band(band: str) -> calibrators.RFCCatalog:
    """Build the full RFC calibrator catalog (no flux cut) for the RFC band matching an observing band.

    Parameters
    ----------
    band : str
        Observing band (e.g. '6cm').

    Returns
    -------
    RFCCatalog
        Catalog to be reused across `_resolve_calibrators` calls for the same band.
    """
    _load_heavy()
    return calibrators.RFCCatalog(min_flux=0.0 * u.Jy, band=calibrators._wavelength_to_rfc_band(band))


def _resolve_calibrators(names: list[str], target: sources.Source, band: str,
                         source_type: sources.SourceType,
                         auto_func: str = 'phasecal',
                         phase_cal_ref: Optional[sources.Source] = None,
                         catalog: Optional[calibrators.RFCCatalog] = None) -> list[sources.Source]:
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
    catalog : RFCCatalog or None
        Pre-built catalog from `_rfc_catalog_for_band(band)`. If None, it is built here.

    Returns
    -------
    list[Source]
        Resolved Source objects with the given source_type assigned. Names that cannot be
        resolved (RFC catalog nor online services) are reported and skipped.
    """
    _load_heavy()
    from astropy.coordinates.name_resolve import NameResolveError
    resolved: list[sources.Source] = []
    rfc_band = calibrators._wavelength_to_rfc_band(band)
    cat = catalog if catalog is not None else _rfc_catalog_for_band(band)

    if not names:
        # Auto-select
        if auto_func == 'phasecal':
            rprint("[yellow]No phase calibrator name given — auto-selecting the best candidate. "
                   "A better source may be found manually via 'planobs phasecals'.[/yellow]")
            src: calibrators.CalibratorSource | None = calibrators.select_phase_calibrator(target, band, catalog=cat)
        else:
            rprint("[yellow]No check source name given — auto-selecting the best candidate. "
                   "A better source may be found manually via 'planobs phasecals'.[/yellow]")
            if phase_cal_ref is None:
                phase_cal_ref = target
            src: calibrators.CalibratorSource | None = calibrators.select_check_source(target, phase_cal_ref, band, catalog=cat)

        if src is not None:
            rprint(f"[green]  → Selected {escape(src.name)} (sep "
                   f"{target.coord.separation(src.coord).deg:.2f}°, "
                   f"unresolved {src.unresolved_flux(rfc_band):.2f} Jy)[/green]")
            resolved.append(sources.Source(name=src.name, coordinates=src.coord,
                                           source_type=source_type, other_names=[src.ivsname]))
        else:
            rprint(f"[bold red]Could not auto-select a {source_type.name.lower()} "
                   f"near {escape(target.name)}.[/bold red]")
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
                    except (ValueError, NameResolveError) as err:
                        rprint(f"[bold red]Source '{escape(sname)}' not found in RFC catalog or external "
                               f"services — skipping ({escape(str(err))}).[/bold red]")

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

    Parameters
    ----------
    band : str
        Observing band, as defined in the catalogs as 'XXcm', with 'XX' being
        the wavelength in cm. See .freqsetups.bands to get a list.
    networks : list of str, optional
        The VLBI network(s) that will participate in the observation.
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
    targets : list[str], optional
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
        Must be > 0 and <= MAX_DURATION_HOURS (96 h).
    datarate : astropy.units.Quantity, optional
        Maximum data rate of the observation (e.g. 4*Gbit/s).
    ontarget : float, optional
        Fraction of the total time of the observation spent on the target source.
        If multiple sources are given (e.g. already specifying scan lengths), this will be ignored.
        Default is 0.7.
    subbands : int, optional
        Number of subbands in which the total bandwidth is split. Default is 4.
    channels : int, optional
        Number of spectral channels in which each subband is divided. Default is 64.
    polarizations : int, optional
        Number of polarizations recorded. It can be 1 (single pol.), 2 (dual pol.), or 4 (full stokes
        recorded: RR, LL, RL, LR). Default is 4.
    inttime : astropy.units.Quantity, optional
        Integration time used in the observations (e.g. time resolution on the correlated data).
        Default is 2 s.
    phasecal_names : list[str], optional
        Phase calibrator source names for the target.
    check_source_names : list[str], optional
        Check source names for the target.
    fringefinder_spec : list[str], optional
        Fringe finder source specification.
    polcal : bool, optional
        Requires polarization calibration for the observation. Default is False.

    Returns
    -------
    VLBIObs
        A VLBI Observation object with all defined parameters.

    Raises
    ------
    ValueError
        If no network/stations are given, the start time is not UTC, the duration is out of
        bounds, a target cannot be resolved, or no valid station is selected.
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

    if duration is not None:
        assert isinstance(duration, u.Quantity)
        duration_hours = duration.to(u.h).value
        if not 0 < duration_hours <= MAX_DURATION_HOURS:
            raise ValueError(f"The observation duration must be > 0 and <= {MAX_DURATION_HOURS:g} hours "
                             f"(got {duration_hours:g} h).")

    src2observe: dict[str, sources.ScanBlock] = {}
    rfc_catalog: Optional[calibrators.RFCCatalog] = None
    if targets is not None:
        from astropy.coordinates.name_resolve import NameResolveError
        for target in targets:
            if (source_catalog is not None) and (target in source_catalog.blocknames):
                src2observe[target] = source_catalog[target]
                continue

            try:
                a_source = sources.Source.source_from_str(target)
            except (ValueError, NameResolveError) as err:
                raise ValueError(
                    f"Source '{target}' not found: not in the provided catalog, not in the "
                    "RFC calibrator catalog, and could not be resolved online "
                    "(SIMBAD/NED/VizieR). Check the source name or add it to your catalog.") from err

            scans_for_block: list[sources.Scan] = []
            if (phasecal_names is not None or check_source_names is not None) and rfc_catalog is None:
                rfc_catalog = _rfc_catalog_for_band(band)

            # Resolve phase calibrator(s)
            if phasecal_names is not None:
                pc_sources = _resolve_calibrators(phasecal_names, a_source, band, sources.SourceType.PHASECAL,
                                                  auto_func='phasecal', catalog=rfc_catalog)
                for pc in pc_sources:
                    scans_for_block.append(sources.Scan(pc, duration=1.5 * u.min))

            # Resolve check source(s)
            if check_source_names is not None:
                # Need a phase cal reference for geometry; use first resolved pc or target
                pc_ref = scans_for_block[0].source if scans_for_block else a_source
                cs_sources = _resolve_calibrators(check_source_names, a_source, band, sources.SourceType.CHECKSOURCE,
                                                  auto_func='check', phase_cal_ref=pc_ref, catalog=rfc_catalog)
                for cs in cs_sources:
                    scans_for_block.append(sources.Scan(cs, duration=1.5 * u.min, every=4))

            # Target scan
            target_dur = freqsetups.phaseref_cycle(band)
            if target_dur is not None and scans_for_block:
                # Subtract phase-cal time from cycle to get target time
                pc_time = sum((s.duration for s in scans_for_block if s.source.type == sources.SourceType.PHASECAL),
                              start=0.0 * u.min)
                target_dur = max(target_dur - pc_time, 1.0 * u.min)
            elif target_dur is None:
                target_dur = 5.0 * u.min
            scans_for_block.append(sources.Scan(a_source, duration=target_dur))
            src2observe[a_source.name] = sources.ScanBlock(scans_for_block)
    elif source_catalog is not None:
        src2observe = source_catalog.blocks
    # No targets and no catalog: valid (e.g. thermal noise for an unspecified source).

    cycle = freqsetups.phaseref_cycle(band)
    default_target_duration = cycle - 1.5*u.min if cycle is not None else 5.0*u.min
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
                        ascan.duration = default_target_duration

    duration_val = duration.to(u.min).value if duration is not None else None
    # Validate networks/stations first: unknown names raise a clean ValueError here.
    selected_stations, excluded_stations = get_stations(band, networks, stations)
    if datarate is None:
        if networks is not None:
            network_station_codes = []
            for n in networks:
                if n in obs._NETWORKS:
                    for s in obs._NETWORKS[n].station_codenames:
                        if s not in network_station_codes:
                            network_station_codes.append(s)
            filtered_stations = obs._STATIONS.filter_antennas(network_station_codes + (stations or []))
        else:
            filtered_stations = obs._STATIONS.filter_antennas(stations or [])
        if not filtered_stations:
            raise ValueError(f"No valid stations provided. Check station codes: {stations}")
        # User-given networks keep their order; guessed ones are sorted by best match first.
        obs_networks = networks or obs.Observation.guess_network(band, filtered_stations)
        known_networks = [n for n in obs_networks if n in obs._NETWORKS]
        if known_networks:
            datarate = obs._NETWORKS[known_networks[0]].max_datarate(band)
    elif isinstance(datarate, int):
        datarate = datarate*u.Mbit/u.s
        rprint("[yellow]Data rate as an int, assumed Mbit/s, but it should have had units[/yellow]")
    elif isinstance(datarate, str):
        raise ValueError("Data rate is a str! ", datarate)

    if start_time is not None and duration_val is not None:
        # 10-min sampling plus the exact end time, strictly increasing and without duplicates.
        times = start_time + np.unique(np.append(np.arange(0, duration_val, 10), duration_val))*u.min
    else:
        times = None
    o = VLBIObs(band, selected_stations, scans=src2observe,
                times=times, duration=duration,
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
    # Handle --get-key-template and --version before full parsing.
    if len(sys.argv) > 1 and sys.argv[1] in ('-V', '--version'):
        print(f"planobs {version('vlbiplanobs')}")
        sys.exit(0)

    if _maybe_handle_get_key_template():
        sys.exit(0)

    # Check if this is legacy mode (no subcommand provided)
    if len(sys.argv) == 1:
        # No arguments - show subcommand help
        parser = argparse.ArgumentParser(description=_MODES_HELP, prog="planobs",
                                         formatter_class=RawTextRichHelpFormatter)
        parser.add_argument('-V', '--version', action='version', version=f"%(prog)s {version('vlbiplanobs')}")
        subparsers = parser.add_subparsers(dest='command', help='Available commands')
        subparsers.add_parser('observe', help='Plan VLBI observations (default mode)')
        subparsers.add_parser('fringefinders', help='Find fringe finder sources')
        subparsers.add_parser('phasecals', help='Find phase calibrator sources')
        subparsers.add_parser('source', help='Get information about a specific source')
        subparsers.add_parser('antenna', aliases=['ant'], help='Get information about antennas')
        subparsers.add_parser('server', help='Start the web server')
        parser.print_help()
        sys.exit(0)
    elif len(sys.argv) > 1 and sys.argv[1] not in ('observe', 'fringefinders', 'phasecals', 'source', 'server', 'antenna', 'ant'):
        # Legacy mode: treat as observation planning
        parser = _parser_with_pager(description=_MODES_HELP, prog="planobs", formatter_class=RawTextRichHelpFormatter)
        parser.add_argument('-V', '--version', action='version', version=f"%(prog)s {version('vlbiplanobs')}")
        add_observation_arguments(parser)
        args = parser.parse_args()
        args.command = 'observe'
        _maybe_setup_logging(getattr(args, 'logging', False))
        handle_observation_command(args)
        return

    parser = _parser_with_pager(description=_MODES_HELP, prog="planobs", formatter_class=RawTextRichHelpFormatter)
    parser.add_argument('-V', '--version', action='version', version=f"%(prog)s {version('vlbiplanobs')}")

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    obs_parser = subparsers.add_parser('observe', add_help=False,
                                       help='Plan VLBI observations (default mode)',
                                       formatter_class=RawTextRichHelpFormatter)
    _add_pager_help(obs_parser)
    add_observation_arguments(obs_parser)

    fringe_parser = subparsers.add_parser('fringefinders', add_help=False,
                                          help='Find fringe finder sources',
                                          formatter_class=RawTextRichHelpFormatter)
    _add_pager_help(fringe_parser)
    add_fringe_finder_arguments(fringe_parser)

    phase_parser = subparsers.add_parser('phasecals', add_help=False,
                                         help='Find phase calibrator sources near a target',
                                         formatter_class=RawTextRichHelpFormatter)
    _add_pager_help(phase_parser)
    add_phase_cal_arguments(phase_parser)

    source_parser = subparsers.add_parser('source', add_help=False,
                                          help='Get information about a specific source',
                                          formatter_class=RawTextRichHelpFormatter)
    _add_pager_help(source_parser)
    add_source_arguments(source_parser)

    server_parser = subparsers.add_parser('server', add_help=False,
                                          help='Start the PlanObs web server',
                                          formatter_class=RawTextRichHelpFormatter)
    _add_pager_help(server_parser)
    add_server_arguments(server_parser)

    antenna_parser = subparsers.add_parser('antenna', aliases=['ant'], add_help=False,
                                           help='Get information about a specific antenna or list antennas by band',
                                           formatter_class=RawTextRichHelpFormatter)
    _add_pager_help(antenna_parser)
    add_antenna_arguments(antenna_parser)

    for subparser in (fringe_parser, phase_parser, source_parser,
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
    """Add grouped arguments for observation planning."""
    main_group = parser.add_argument_group('Main parameters')
    main_group.add_argument('-b', '--band', type=str, help="Observing band, as defined "
                            "in the catalogs as 'XXcm', with 'XX' being\nthe wavelegnth in cm. "
                            "[green]See '--list-bands' to get a list.[/green]")
    main_group.add_argument('-n', '--network', type=str, nargs='+',
                            help="The VLBI network(s) that will participate in\nthe observation. "
                            "It will take the default stations in each network.\nIf 'stations' "
                            "is provided, then it will take both the default stations\nplus "
                            "the ones given in stations. [green]See '--list-networks' to get a list.[/green]")
    main_group.add_argument('-s', '--stations', type=str, nargs='+', help="List "
                            "of the antennas that will participate in the\nobservation. "
                            "You can use either antenna codenames or the standard name,\n"
                            "as given in the catalogs. [green]See '--list-antennas' to get a list.[/green]")
    main_group.add_argument('-e', '--epoch', type=str, default=None, dest='epoch',
                            help="Start of the observation, with the format 'YYYY-MM-DD HH:MM' "
                            "in UTC.")
    add_deprecated_alias(main_group, '-t1', '--starttime', dest='epoch', new_option='-e/--epoch', type=str)
    main_group.add_argument('-d', '--duration', type=float, default=None,
                            help="Total duration of the observation, in hours.")

    source_group = parser.add_argument_group('Source-related options')
    source_group.add_argument('-t', '--target', type=str, default=None, nargs='+', dest='targets', metavar='TARGET',
                              help="Source(s) to be observed. Each entry can be:\n"
                              "  a) A source name (looked up in SIMBAD/NED/VizieR/RFC).\n"
                              "  b) Coordinates: 'hh:mm:ss dd:mm:ss' or 'XXhXXmXXs XXdXXmXXs'.\n"
                              "  c) 'name/coordinates' to provide both (coordinates override lookup).\n"
                              "Or if '--source-catalog' is defined, selects the block(s) in that file.\n"
                              "Multiple sources can be provided.")
    add_deprecated_alias(source_group, '--targets', dest='targets', new_option='-t/--target', type=str, nargs='+')
    source_group.add_argument('-sc', '--source-catalog', '--sc', type=str, default=None,
                              help="Input file containing the personal source catalog.\n"
                              "If provided, then '--target' will select the block(s) "
                              "defined in\nthis file, ignoring the rest.")
    source_group.add_argument('--fringefinders', default=['2'], type=str, nargs='+',
                                help="Defines the fringe finder source(s) to be scheduled "
                                "in the observation.\nIt can be either a list of source names "
                                "(as long as they\nappear in AstroGeo), "
                                "'name/coordinates' to provide both,\n"
                                "or a single number, meaning how many scans should go on\nfringe "
                                "finders, and it will automatically select the most suitable sources.")
    source_group.add_argument('--polcal', action="store_true", default=False,
                              help="Requires polarization calibration for the observation.")
    source_group.add_argument('--phasecal', default=None, type=str, nargs='*',
                              help="Phase calibrator source(s) for the target. If no names are given\n"
                              "(just --phasecal), the best candidate is picked automatically.\n"
                              "One or more source names or 'name/coordinates' can be provided.")
    source_group.add_argument('--check-source', default=None, type=str, nargs='*',
                              help="Check source(s) for the target. If no names are given\n"
                              "(just --check-source), the best candidate is picked automatically.\n"
                              "One or more source names or 'name/coordinates' can be provided.")
    source_group.add_argument('--pulsar', default=None, type=str,
                              help="Sets to schedule at least a scan on a pulsar source. "
                              "If a number,\nit will select a pulsar from the personal "
                              "input source file (must be\nprovided!). If a name, 'name/coordinates', "
                              "or coordinates, it will\nresolve accordingly.")

    antenna_group = parser.add_argument_group('Antenna-related options')
    antenna_group.add_argument('--station-catalog', type=str, default=None,
                               help="Input file containing the personal station catalog.\n"
                               "If provided, then the default catalog will not be read.")
    antenna_group.add_argument('--list-antennas', action="store_true", default=False,
                               help="Prints the list of all antennas defined in PlanObs.")
    antenna_group.add_argument('--list-networks', action="store_true", default=False,
                               help="Prints the list of all VLBI networks defined in PlanObs.")

    obs_group = parser.add_argument_group('Observation configuration')
    obs_group.add_argument('--list-bands', action="store_true", default=False,
                           help="Writes the list of all observing bands defined in PlanObs.")
    obs_group.add_argument('--data-rate', type=float, default=None,
                           help="Maximum data rate of the observation, in Mb/s.")
    obs_group.add_argument('--debug', action="store_true", default=False,
                           help="If set, shows some debuging messages.")
    obs_group.add_argument('--nme', action="store_true", default=False,
                           help="Network Monitoring Experiment mode (no targets needed; requires\n"
                           "'-e' and '-d'). The full time is covered with ~15-min fringe-finder\n"
                           "scans visible by all antennas, with ftp fringe-test grabs every 30 min\n"
                           "(every 15 min if the duration is <= 2.5 h). Without '--sched' it lists\n"
                           "the visible fringe finders and the proposed scans; with '--sched' it\n"
                           "writes the NME .key file. Explicit '--fringefinders' names are used\n"
                           "as the only candidates.")

    output_group = parser.add_argument_group('Output options')
    output_group.add_argument('--sched', default=None, type=str,
                              help="Produces a (SCHED) .key schedule file for "
                              "the observation with the\ngiven name.")
    output_group.add_argument('--setup', default=None, type=str,
                              help="Frequency setup to write in the 'setup = ...' line of the .key\n"
                              "file produced by --sched. If not given, PlanObs guesses it from the\n"
                              "observation setup.")
    output_group.add_argument('--template', default=None, type=str,
                              help="SCHED .key template file. If not given, the default template\n"
                              "distributed with PlanObs is used.")
    output_group.add_argument('-o', '--output', type=str, default=None, metavar='FILENAME',
                              help="Write all observation inputs and results to .pdf, .txt, .md, or .json.\n"
                                   "The filename extension selects the format (case-insensitive).")
    add_logging_argument(output_group)


def add_fringe_finder_arguments(parser):
    """Add arguments for fringe finder search (also used by calibrators.main_fringe)."""
    parser.add_argument('-n', '--network', type=str, nargs='+',
                        help="The VLBI network(s) that will participate in\nthe observation. "
                        "It will take the default stations in each network.\nIf 'stations' "
                        "is provided, then it will take both the default stations\nplus "
                        "the ones given in stations. [green]See '--list-networks' to get a list.[/green]")
    parser.add_argument('-s', '--stations', type=str, nargs='+',
                        help="List of antenna codenames or names that will participate in the observation.")
    # Not required=True at the argparse level so that the deprecated aliases also satisfy it; enforced in
    # handle_fringe_finder_command.
    parser.add_argument('-e', '--epoch', type=str, default=None, dest='epoch',
                        help="Start of the observation in format 'YYYY-MM-DD HH:MM' (UTC). Required.")
    add_deprecated_alias(parser, '-t1', '--starttime', dest='epoch', new_option='-e/--epoch', type=str)
    add_removed_option(parser, '-t', new_option='-e/--epoch', reason="'-t' means '--target' in the other commands")
    parser.add_argument('-d', '--duration', type=float, required=True, help="Duration of the observation in hours.")
    # Defaults are None here and resolved in handle_fringe_finder_command from the calibrators.FRINGE_DEFAULT_*
    # constants (single source of truth); calibrators is heavy and is only imported lazily (see _load_heavy).
    # The '(default: X)' help texts are checked against those constants in tests/test_fixes_cli.py.
    parser.add_argument('--min-flux', type=float, default=None,
                        help="Minimum unresolved flux threshold in Jy (default: 0.5).")
    parser.add_argument('--min-elevation', type=float, default=None,
                        help="Minimum elevation in degrees (default: 20).")
    parser.add_argument('-l', '--max-lines', type=int, default=None,
                        help="Maximum number of sources to return (default: 20).")
    parser.add_argument('--require-all', action='store_true', default=False,
                        help="Require source to be visible by ALL stations (default: False).")
    parser.add_argument('-b', '--band', type=str, default=None,
                        help="Observing band for flux display (e.g., '18cm', '6cm'). If not provided, shows flux for all available bands.")
    parser.add_argument('--station-catalog', type=str, default=None, help="Path to custom station catalog file.")
    parser.add_argument('--json', action='store_true', default=False,
                        help="Output results in JSON format instead of a table.")


def add_phase_cal_arguments(parser):
    """Add arguments for phase calibrator search (also used by calibrators.main_phasecal)."""
    parser.add_argument('target', type=str, nargs='?', default=None,
                        help="Target source name (J2000 or IVS name from RFC catalog;\n"
                        "or a block/source name from '--source-catalog' if provided).")
    parser.add_argument('-t', '--target', type=str, default=None, dest='target_option', metavar='TARGET',
                        help="Target source, as an alternative to the positional TARGET argument.")
    parser.add_argument('-sc', '--source-catalog', '--sc', type=str, default=None,
                        help="Input file containing the personal source catalog.\n"
                        "If provided, then '--target' will first be looked up in\nthis file "
                        "(block or source name), as in the observe mode.")
    # Defaults are None here and resolved in handle_phase_cal_command from the calibrators.PHASECAL_DEFAULT_*
    # constants (see add_fringe_finder_arguments for the rationale).
    parser.add_argument('--max-separation', type=float, default=None,
                        help="Maximum angular separation in degrees (default: 5.0).")
    parser.add_argument('--min-flux', type=float, default=None,
                        help="Minimum unresolved flux threshold in Jy (default: 0.1).")
    parser.add_argument('-l', '--max-lines', type=int, default=None, dest='max_lines',
                        help="Maximum number of sources to return (default: all).")
    add_deprecated_alias(parser, '--n-sources', dest='max_lines', new_option='-l/--max-lines', type=int)
    add_removed_option(parser, '-n', new_option='-l/--max-lines', reason="'-n' means '--network' in the other commands")
    parser.add_argument('-b', '--band', type=str, default=None,
                        help="Observing band for flux display (e.g., '18cm', '6cm'). If not provided, shows flux for all available bands.")
    parser.add_argument('--rfc-catalog', type=str, default=None, dest='rfc_catalog',
                        help="Path to custom RFC catalog file.")
    add_deprecated_alias(parser, '--catalog-file', dest='rfc_catalog', new_option='--rfc-catalog', type=str)
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


def add_logging_argument(parser_or_group):
    """Add the optional ``--logging`` flag that enables writing a log file.

    No log file is created by default. When the flag is given without a value,
    a default path is used ('/var/log/planobs.log' if writable, otherwise
    '~/log-planobs.log'). An explicit path may also be provided.

    Parameters
    ----------
    parser_or_group : argparse.ArgumentParser or argparse._ArgumentGroup
        Parser or argument group to which the flag should be added.
    """
    parser_or_group.add_argument('--logging', nargs='?', const=True, default=False, metavar='LOGFILE',
                                 help="Enable logging to a file. Optionally provide a path;\n"
                                 "otherwise '/var/log/planobs.log' (if writable) or\n"
                                 "'~/log-planobs.log' is used. Disabled by default.")


def _sched_paths(sched: str) -> tuple[str, str]:
    """Derive the .key output filename and the SCHED experiment code from the --sched argument.

    Parameters
    ----------
    sched : str
        Value given to --sched (with or without the '.key' extension, may include directories).

    Returns
    -------
    tuple[str, str]
        (key_filename, experiment_code). The code is the file stem, uppercased.

    Raises
    ------
    ValueError
        If the file stem is not a valid experiment code (only letters, digits and '_').
    """
    key_filename = sched if sched.endswith('.key') else f"{sched}.key"
    experiment_code = Path(key_filename).stem.upper()
    if not re.fullmatch(_EXPERIMENT_CODE_PATTERN, experiment_code):
        raise ValueError(f"Invalid experiment code '{experiment_code}' derived from --sched '{sched}': "
                         "the file name may only contain letters, digits and '_' (e.g. 'EM179A.key').")
    return key_filename, experiment_code


def _write_key_file(key_filename: str, key_content: str, label: str = 'Schedule') -> None:
    """Write a SCHED .key file and report it to the user.

    Parameters
    ----------
    key_filename : str
        Destination path (overwritten if it exists).
    key_content : str
        Full .key file content.
    label : str
        Human-readable label used in the confirmation message (e.g. 'Schedule', 'NME schedule').
    """
    with open(key_filename, 'w') as f:
        f.write(key_content)
    rprint(f"[green]{label} file written to: {escape(key_filename)}[/green]")


def handle_observation_command(args):
    """Handle the observation planning command (also dispatches the --nme mode).

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments from `add_observation_arguments`. Exits the process with status 1 on
        invalid input or failures.
    """
    _load_heavy()
    t0 = dt.now() if args.debug else None

    if args.list_networks:
        rprint("[bold]Available VLBI networks:[/bold]")
        for network_name, network in obs._NETWORKS.items():
            rprint(f"[bold]{escape(network_name)}:[/bold] {escape(network.name)}")
            rprint(f"  [dim]Default antennas: {escape(', '.join(network.station_codenames))}[/dim]")
            rprint(f"  [dim]Observes at: {', '.join(network.observing_bands)}[/dim]")

    if args.list_antennas:
        rprint("\n[bold]All available antennas:[/bold]")
        for ant in obs._STATIONS:
            rprint(f"     {escape(ant.name)} ({escape(ant.codename)}):  {ant.diameter} in {escape(ant.country)}")
            rprint(f"      [dim]Observes at {', '.join(ant.bands)}[/dim]")

    if args.list_bands:
        rprint("\n[bold]Available observing bands:[/bold]")
        for aband in obs.freqsetups.bands:
            rprint(f"[bold]{aband}[/bold] [dim]({obs.freqsetups.bands[aband]})[/dim]")
            rprint("[dim]  Observable with [/dim]", end='')
            rprint("[dim]" +
                   escape(', '.join([nn for nn, n in obs._NETWORKS.items() if aband in n.observing_bands])) +
                   "[/dim]")

    if args.list_antennas or args.list_bands or args.list_networks:
        sys.exit(0)

    if args.band is None:
        rprint("\n\n[bold red]The observing band (-b/--band) and either '--network' and/or "
               "'--stations' are mandatory.[/bold red]")
        sys.exit(1)

    if args.band not in obs.freqsetups.bands:
        rprint(f"[bold red]The provided band ({escape(args.band)}) is not available"
               "[/bold red]\n[dim]These are the available bands: "
               f"{', '.join(obs.freqsetups.bands)}.[/dim]")
        sys.exit(1)

    if args.network is None and args.stations is None:
        rprint("[bold red]You need to provide at least a VLBI network "
               "or a list of antennas that will participate in the observation.[/bold red]")
        sys.exit(1)

    if args.duration is None and args.epoch is not None:
        rprint("[bold red]If you provide a start time, you also need to provide a duration for"
               " the observation.[/bold red]")
        sys.exit(1)

    if args.data_rate is not None and args.data_rate <= 0:
        rprint(f"[bold red]The data rate must be a positive number in Mb/s (got {args.data_rate:g}).[/bold red]")
        sys.exit(1)

    if args.sched is not None:
        try:
            _sched_paths(args.sched)
        except ValueError as error:
            rprint(f"[bold red]Error: {escape(str(error))}[/bold red]")
            sys.exit(1)

    if getattr(args, 'nme', False):
        handle_nme_command(args)
        return

    if getattr(args, 'setup', None) is not None and args.sched is None:
        rprint("[bold yellow]--setup is only used when producing a schedule file (--sched). "
               "Ignoring it.[/bold yellow]")

    output_filename = getattr(args, 'output', None)
    if output_filename is not None:
        from vlbiplanobs.report import validate_report_filename
        try:
            validate_report_filename(output_filename)
        except ValueError as error:
            rprint(f"[bold red]Error: {escape(str(error))}[/bold red]")
            sys.exit(1)

    # Resolve phasecal / check-source arguments (None = not requested, [] = auto-select)
    phasecal_arg = getattr(args, 'phasecal', None)
    check_source_arg = getattr(args, 'check_source', None)
    fringefinder_arg = getattr(args, 'fringefinders', ['2'])
    polcal_arg = getattr(args, 'polcal', False)

    try:
        o = main(band=args.band, networks=args.network, stations=args.stations,
                 src_catalog=args.source_catalog, station_catalog=args.station_catalog,
                 targets=args.targets, start_time=Time(args.epoch, scale='utc') if args.epoch else None,
                 duration=args.duration*u.hour if args.duration is not None else None,
                 datarate=args.data_rate*u.Mbit/u.s if args.data_rate is not None else None,
                 phasecal_names=phasecal_arg, check_source_names=check_source_arg,
                 fringefinder_spec=fringefinder_arg, polcal=polcal_arg)
    except ValueError as e:
        rprint(f"[bold red]Error: {escape(str(e))}[/bold red]")
        sys.exit(1)

    o.summary(gui=False, tui=True)
    if args.targets is not None or args.source_catalog is not None:
        o.plot_visibility(gui=False, tui=True)

    if output_filename is not None:
        from vlbiplanobs.report import write_observation_report
        try:
            output_path = write_observation_report(o, output_filename)
        except Exception as error:
            rprint(f"[bold red]Could not write observation report: {escape(str(error))}[/bold red]")
            sys.exit(1)
        rprint(f"[green]Observation report written to: {escape(str(output_path))}[/green]")

    if args.sched is not None:
        key_filename, experiment_code = _sched_paths(args.sched)
        from vlbiplanobs.scheduler import ObservationScheduler
        scheduler = ObservationScheduler(o, fringefinder_spec=fringefinder_arg, polcal=polcal_arg)
        scheduler.schedule()
        try:
            key_content = scheduler.generate_key_file(experiment_code=experiment_code,
                                                      setup_file=getattr(args, 'setup', None),
                                                      template_path=getattr(args, 'template', None))
        except ValueError as error:
            rprint(f"[bold red]Could not generate schedule file: {escape(str(error))}[/bold red]")
            sys.exit(1)
        _write_key_file(key_filename, key_content, label='Schedule')

    if args.debug:
        print(f"Execution time: {(dt.now() - t0).total_seconds()} s")


def _print_nme_candidates(scans, sources, counts, times, n_stations: int, band: str, max_lines: int = 20) -> None:
    """Print the fringe finders visible by all antennas at some point of the NME, best coverage first.

    Parameters
    ----------
    scans : list[nme.NMEScan]
        Planned scans (used to mark the sources selected for the schedule).
    sources : list[Source]
        Candidate sources, columns of ``counts``.
    counts : np.ndarray
        Number of stations observing each source at each time, shape (n_times, n_sources).
    times : Time
        Sampling times of ``counts``.
    n_stations : int
        Number of participating stations.
    band : str
        Observing band, used for the flux column.
    max_lines : int
        Maximum number of candidates to print.
    """
    from vlbiplanobs import nme
    full = counts >= n_stations
    coverage = full.mean(axis=0)
    order = [i for i in np.argsort(-coverage, kind='stable') if coverage[i] > 0]
    selected = {s.source.name for s in scans}
    rprint(f"\n[bold green]{len(order)} fringe finder candidates are visible by all {n_stations} "
           f"antennas at some point of the observation:[/bold green]")
    table = Table(show_header=True, header_style="bold", box=box.SIMPLE)
    for col, justify in (("Name", "left"), ("IVS Name", "left"), ("Unresolved (Jy)", "right"),
                         ("All-antenna time", "right"), ("UTC windows (all antennas)", "left"), ("In plan", "center")):
        table.add_column(col, justify=justify)
    for i in order[:max_lines]:
        src = sources[i]
        windows = ', '.join(f"{a.datetime.strftime('%H:%M')}-{b.datetime.strftime('%H:%M')}"
                            for a, b in nme.visible_windows(full[:, i], times))
        flux = nme.source_flux(src, band)
        table.add_row(escape(src.name), escape(getattr(src, 'ivsname', '') or ''), f"{flux:.2f}" if flux > 0 else "N/A",
                      f"{100 * coverage[i]:.0f}%", windows, "*" if src.name in selected else "")
    rprint(table)
    if len(order) > max_lines:
        rprint(f"[dim]... and {len(order) - max_lines} more sources.[/dim]")


def _print_nme_plan(scans, start: Time, n_stations: int) -> None:
    """Print the proposed NME scan list (UTC times, source, visible antennas, grab time)."""
    rprint("\n[bold green]Proposed NME scans:[/bold green]")
    table = Table(show_header=True, header_style="bold", box=box.SIMPLE)
    for col in ("Scan", "Start (UTC)", "Stop (UTC)", "Source", "Antennas", "ftp grab (UTC)"):
        table.add_column(col)
    for k, scan in enumerate(scans, start=1):
        grab = (start + scan.grab_s*u.s).datetime.strftime('%H:%M:%S') if scan.grab_s is not None else ""
        ants = f"{scan.n_visible}/{n_stations}"
        table.add_row(str(k), (start + scan.rec_start_s*u.s).datetime.strftime('%H:%M:%S'),
                      (start + scan.stop_s*u.s).datetime.strftime('%H:%M:%S'), escape(scan.source.name),
                      ants if scan.n_visible == n_stations else f"[bold red]{ants}[/bold red]", grab)
    rprint(table)


def handle_nme_command(args):
    """Handle the Network Monitoring Experiment (--nme) mode of the observe command.

    Plans fringe-finder scans covering the full observation with periodic ftp grabs.
    Without ``--sched``, prints the visible fringe finders and the proposed scans.
    With ``--sched``, also writes the NME SCHED .key file.
    """
    from vlbiplanobs import nme
    from vlbiplanobs.scheduler import format_setup_line, guess_setup_file
    if args.epoch is None or args.duration is None:
        rprint("[bold red]The NME mode requires both the start epoch (-e/--epoch) and the duration (-d).[/bold red]")
        sys.exit(1)
    if args.targets is not None or args.source_catalog is not None:
        rprint("[bold yellow]Targets/source catalogs are ignored in NME mode.[/bold yellow]")

    ff_arg = getattr(args, 'fringefinders', ['2'])
    ff_names = None if (len(ff_arg) == 1 and ff_arg[0].isdigit()) else ff_arg
    try:
        o = main(band=args.band, networks=args.network, stations=args.stations,
                 station_catalog=args.station_catalog, start_time=Time(args.epoch, scale='utc'),
                 duration=args.duration*u.hour,
                 datarate=args.data_rate*u.Mbit/u.s if args.data_rate is not None else None)
        fringefinders = nme.resolve_fringe_finders(ff_names) if ff_names else None
        start = Time(args.epoch, scale='utc')
        scans, sources_ff, counts, times = nme.plan_nme(o.stations, start, args.duration*u.hour,
                                                        args.band, fringefinders=fringefinders)
    except ValueError as e:
        rprint(f"[bold red]Error: {escape(str(e))}[/bold red]")
        sys.exit(1)

    n_stations = len(o.stations)
    rprint(f"[bold]NME at {args.band} with {n_stations} antennas: "
           f"{', '.join(s.codename for s in o.stations)}[/bold]")
    if o._excluded_stations:
        rprint(f"[yellow]Antennas dropped (cannot observe at {args.band}): "
               f"{', '.join(o._excluded_stations)}[/yellow]")
    _print_nme_candidates(scans, sources_ff, counts, times, n_stations, args.band)
    _print_nme_plan(scans, start, n_stations)
    if any(s.n_visible < n_stations for s in scans):
        rprint("[bold yellow]Some scans cannot be observed by all antennas (see the table above).[/bold yellow]")

    if args.sched is None:
        return

    key_filename, experiment_code = _sched_paths(args.sched)
    setup_file = getattr(args, 'setup', None) or guess_setup_file(o)
    datarate = int(o.datarate.to(u.Mbit/u.s).value) if o.datarate is not None else None
    try:
        key_content = nme.generate_nme_key_file(scans, o.stations, start, args.band, experiment_code,
                                                format_setup_line(setup_file), datarate_mbps=datarate,
                                                template_path=getattr(args, 'template', None))
    except (ValueError, OSError) as error:
        rprint(f"[bold red]Could not generate the NME schedule file: {escape(str(error))}[/bold red]")
        sys.exit(1)
    _write_key_file(key_filename, key_content, label='NME schedule')


def handle_fringe_finder_command(args):
    """Handles 'planobs fringefinders': validates the parsed args and runs calibrators.run_fringe_finders.

    Inputs
        args : argparse.Namespace — parsed arguments from add_fringe_finder_arguments.

    Exits
        With the exit code returned by calibrators.run_fringe_finders (1 if no epoch or no network/stations given).
    """
    if args.epoch is None:
        rprint("[bold red]The start of the observation (-e/--epoch 'YYYY-MM-DD HH:MM') is required.[/bold red]")
        sys.exit(1)
    _load_heavy()
    if args.network is None and args.stations is None:
        rprint("[bold red]You need to provide at least a VLBI network "
               "or a list of antennas that will participate in the observation.[/bold red]")
        sys.exit(1)

    min_flux = calibrators.FRINGE_DEFAULT_MIN_FLUX_JY if args.min_flux is None else args.min_flux
    min_elevation = calibrators.FRINGE_DEFAULT_MIN_ELEVATION_DEG if args.min_elevation is None else args.min_elevation
    max_lines = calibrators.FRINGE_DEFAULT_MAX_LINES if args.max_lines is None else args.max_lines
    sys.exit(calibrators.run_fringe_finders(starttime=args.epoch, duration=args.duration, networks=args.network,
                                            stations=args.stations, min_flux=min_flux, min_elevation=min_elevation,
                                            max_lines=max_lines, require_all=args.require_all, band=args.band,
                                            station_catalog=args.station_catalog, as_json=args.json))


def handle_phase_cal_command(args):
    """Handles 'planobs phasecals': resolves the target (positional or '-t/--target') and runs
    calibrators.run_phasecals.

    Inputs
        args : argparse.Namespace — parsed arguments from add_phase_cal_arguments.

    Exits
        With the exit code returned by calibrators.run_phasecals (1 if the target is missing or ambiguous).
    """
    target = args.target if args.target is not None else args.target_option
    if target is None:
        rprint("[bold red]A target source is required: 'planobs phasecals TARGET \\[options]'.[/bold red]")
        sys.exit(1)
    if args.target is not None and args.target_option is not None and args.target != args.target_option:
        rprint("[bold red]Two different targets were given (positional and '-t/--target'); "
               "provide only one.[/bold red]")
        sys.exit(1)

    _load_heavy()
    max_separation = calibrators.PHASECAL_DEFAULT_MAX_SEPARATION_DEG if args.max_separation is None \
        else args.max_separation
    min_flux = calibrators.PHASECAL_DEFAULT_MIN_FLUX_JY if args.min_flux is None else args.min_flux
    sys.exit(calibrators.run_phasecals(target=target, max_separation=max_separation, min_flux=min_flux,
                                       n_sources=args.max_lines, band=args.band, catalog_file=args.rfc_catalog,
                                       source_catalog=args.source_catalog, as_json=args.json))


def _format_gst_ranges(gst_pairs: list[tuple]) -> str:
    """Format a list of GST (Longitude) start/end pairs as 'HH:MM-HH:MM, ...' string.

    Parameters
    ----------
    gst_pairs : list[tuple]
        List of (start, end) Time pairs.

    Returns
    -------
    str
        Formatted string like 'HH:MM-HH:MM, HH:MM-HH:MM, ...'.
    """
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

    Parameters
    ----------
    args : tuple
        (net_key, band, source_name, ra_deg, dec_deg[, return_gst]).

    Returns
    -------
    tuple
        (net_key, band, is_observable, gst_ranges_str or None).
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
        # GST ranges are the same intervals as the UTC ones, so a single call answers both.
        observable_ranges = observation.when_is_observable(min_stations=3, return_gst=return_gst).get(src.name, [])
        is_obs = bool(observable_ranges)
        if return_gst:
            return net_key, band, is_obs, _format_gst_ranges(observable_ranges) if is_obs else ""
        return net_key, band, is_obs, None
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

    Parameters
    ----------
    source : sources.Source
        The source to check observability for.
    show_gst : bool, optional
        If True, appends a GST column with time ranges where >3 antennas observe.
        Default is False.
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
            rprint(f"[bold]{escape(calibrator.name)}[/bold] (also known as {escape(str(calibrator.ivsname))})")
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
                rprint(f"[bold]Name:[/bold] {escape(source_obj.name)}")
                rprint(f"[bold]Coordinates:[/bold] {source_obj.coord.to_string('hmsdms')}")
            except Exception as e:
                rprint(f"[bold red]Source '{escape(args.source_name)}' not recognized[/bold red]")
                rprint(f"[red]Error: Not found in the RFC Catalog and {escape(str(e))} [/red]")
                sys.exit(1)

        # Show observability table for all sources
        if source_obj and not args.no_networks:
            rprint("\n[bold green]Observable by (bands in cm)[/bold green]")
            _show_observability_table(source_obj, show_gst=args.gst)
    except Exception as e:
        rprint(f"[bold red]Error:[/bold red] {escape(str(e))}")
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

    Parameters
    ----------
    bands_to_show : list[str]
        List of band strings to display (e.g. ['6cm', '18cm']).
    ant : Station
        Station object with sefd(band) method.
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

    Parameters
    ----------
    band : str
        Band string to normalize.

    Returns
    -------
    str
        Normalized band string with 'cm' suffix.
    """
    band = band.strip()
    if band.replace('.', '').isdigit():
        band = f"{band}cm"
    return band


def _find_antenna(name: str) -> Optional[object]:
    """Search obs._STATIONS for an antenna matching name, fullname, or codename (case-insensitive).

    Parameters
    ----------
    name : str
        Antenna name, fullname, or codename to search for.

    Returns
    -------
    Station or None
        The matching Station object or None if not found.
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

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    """
    _load_heavy()
    band = _normalize_band(args.band) if args.band else None
    if band is not None and band not in obs.freqsetups.bands:
        rprint(f"[bold red]Band '{escape(band)}' is not recognized.[/bold red] "
               f"Available bands: {', '.join(obs.freqsetups.bands)}")
        sys.exit(1)

    if args.antenna_name is None:
        if band is None:
            # Same as --list-antennas
            rprint("\n[bold]All available antennas:[/bold]")
            for ant in obs._STATIONS:
                rprint(f"     [bold]{escape(ant.name)}[/bold] ([cyan]{escape(ant.codename)}[/cyan]):  "
                       f"{escape(str(ant.diameter))} in {escape(str(ant.country))}")
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
        rprint(f"[bold red]Antenna '{escape(args.antenna_name)}' not found.[/bold red] "
               "Run [bold]planobs --list-antennas[/bold] to see the available antennas.")
        sys.exit(1)

    rprint(f"\n[bold underline]{escape(ant.fullname)}[/bold underline]")
    rprint(f"  [bold]Short name:[/bold]  {escape(ant.name)}")
    rprint(f"  [bold]Codename:[/bold]    [cyan]{escape(ant.codename)}[/cyan]")
    rprint(f"  [bold]Diameter:[/bold]    {escape(str(ant.diameter))}")
    rprint(f"  [bold]Country:[/bold]     {escape(str(ant.country))}")

    if band is not None:
        if not ant.has_band(band):
            rprint(f"\n[bold red]{escape(ant.name)} does not observe at {band}.[/bold red] "
                   f"It observes at: {', '.join(ant.bands)}")
            sys.exit(1)
        bands_to_show = [band]
    else:
        bands_to_show = sorted(ant.bands, key=lambda b: float(b.replace('cm', '')))

    # Print SEFD table (horizontal layout, wraps at terminal width)
    rprint("")
    _render_horizontal_band_table(bands_to_show, ant)


if __name__ == '__main__':
    cli()
