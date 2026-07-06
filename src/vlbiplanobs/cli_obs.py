"""Heavy CLI observation objects (VLBIObs) and helpers.

This module holds the parts of the CLI that require the heavy dependencies
(numpy, astropy, the vlbiplanobs computation modules). It is imported lazily
from `vlbiplanobs.cli` so that `planobs -h`/`-V`/no-arguments respond immediately.
"""
import sys
import numpy as np
from astropy import units as u
from astropy.time import Time
from rich import print as rprint
from vlbiplanobs import observation as obs
from vlbiplanobs import sources


def optimal_units(value: u.Quantity, units: list[u.Unit]):
    """Given a value (with some units), returns the unit choice from all
    `units` possibilities that better suits the value.

    It is meant for the following use: given 0.02*u.Jy and
    units = [u.kJy, u.Jy, u.mJy, u.uJy], it will return 20*u.mJy.

    Parameters
    ----------
    value : astropy.units.Quantity
        The value to convert.
    units : list[astropy.units.Unit]
        Candidate units, in decreasing scale, all compatible with `value`'s units.

    Returns
    -------
    astropy.units.Quantity
        `value` converted to whichever unit in `units` best fits (a magnitude
        between 0.8 and 800), falling back to the largest or smallest unit if
        `value` is too high or too low for any of them.
    """
    for a_unit in units:
        if 0.8 < value.to(a_unit).value <= 800:
            return value.to(a_unit)

    # Value too high or too low
    if value.to(units[0]).value > 1:
        return value.to(units[0])

    return value.to(units[-1])


class VLBIObs(obs.Observation):
    """CLI-facing extension of `observation.Observation` that adds terminal/GUI
    summary printing, per-source elevation plots (via plotext), and tracking of
    stations that were requested but excluded from the observation.
    """

    def __init__(self, *args, **kwargs) -> None:
        """Initializes a VLBIObs, forwarding all arguments to `observation.Observation`.

        Parameters
        ----------
        *args
            Forwarded unchanged to `observation.Observation.__init__`.
        **kwargs
            Forwarded unchanged to `observation.Observation.__init__`.
        """
        super().__init__(*args, **kwargs)
        self._excluded_stations: dict[str, str] = {}

    def summary(self, gui: bool = True, tui: bool = True):
        """Prints the observation summary. Uses the GUI variant when `gui` is set,
        otherwise the terminal (TUI) variant when `tui` is set; does nothing if both are False
        (the CLI warns about that combination before calling).

        Parameters
        ----------
        gui : bool, optional
            If True, calls the GUI summary variant. Default is True.
        tui : bool, optional
            If True (and `gui` is False), calls the terminal (TUI) summary variant. Default is True.

        Returns
        -------
        None
        """
        if gui:
            return self._summary_gui()

        if tui:
            return self._summary_tui()

        return None

    def _summary_gui(self):
        """GUI variant of the observation summary.

        TODO: for now it just calls the TUI variant (`_summary_tui`); no dedicated
        GUI rendering exists yet.
        """
        self._summary_tui()

    def per_source_elevations(self) -> dict[str, dict[str, np.ndarray]]:
        """Returns per-source elevation data across all stations in degrees.

        Returns
        -------
        dict[str, dict[str, np.ndarray]]
            Elevation values in degrees for each source name and station codename.
        """
        elevations = self.elevations()
        result = {}
        for source_name, station_elevs in elevations.items():
            result[source_name] = {}
            for station_codename, elev_values in station_elevs.items():
                # Convert to degrees and handle astropy Quantities
                if hasattr(elev_values, 'to'):
                    result[source_name][station_codename] = elev_values.to(u.deg).value
                else:
                    result[source_name][station_codename] = np.asarray(elev_values)
        return result

    def _plot_source_elevation_terminal(self, src_name: str, src_vis: dict, src_elev: dict, time_labels: list[str]) -> None:
        """Render a per-antenna elevation plot for one source using plotext.

        Each visible time step is plotted as a colored scatter point.
        Color encodes elevation: red (<10°), yellow (10-20°), green (20-40°),
        cyan (40-60°), blue (>60°). Invisible steps are left blank.

        Parameters
        ----------
        src_name : str
            Source name used as the plot title.
        src_vis : dict[str, np.ndarray]
            Boolean visibility array per antenna codename.
        src_elev : dict[str, np.ndarray]
            Elevation in degrees per antenna codename.
        time_labels : list[str]
            Time label string for every time step (x-axis tick labels).
        """
        import plotext as pltx
        import shutil

        antenna_names = list(src_vis.keys())
        n_times = len(next(iter(src_vis.values())))
        n_ants = len(antenna_names)
        if n_times == 0 or n_ants == 0:
            return

        pltx.clf()
        pltx.theme('clear')

        # Elevation color bands — same scheme as GUI.
        # 256-colour index 226 = bright yellow (avoids 'yellow' rendering as white).
        bands = [
            (0,  10, 'red',  '< 10°'),
            (10, 20, 226,    '10–20°'),
            (20, 40, 'green','20–40°'),
            (40, 60, 'cyan', '40–60°'),
            (60, 90, 'blue', '> 60°'),
        ]

        # Legend width = "▪▪ " prefix (3) + longest label + 1 padding.
        # Extend xlim left by this amount so the legend occupies negative x space
        # and never overlaps the data.  Widen the canvas by the same amount so
        # data density stays at ~1 column per time step.
        legend_offset = max(len(label) for *_, label in bands) + 4
        y_margin = max(len(a) for a in antenna_names) + 4
        term_width = shutil.get_terminal_size((120, 24)).columns
        canvas_width = min(term_width, n_times + y_margin + legend_offset)
        pltx.plot_size(canvas_width, n_ants + 6)

        for el_min, el_max, color, label in bands:
            xs: list[int] = []
            ys: list[int] = []
            for y_idx, ant in enumerate(antenna_names):
                visibility = src_vis.get(ant, np.zeros(n_times, dtype=bool))
                elevations = src_elev.get(ant, np.zeros(n_times))
                for x_idx, (is_vis, elev) in enumerate(zip(visibility, elevations)):
                    if is_vis and el_min <= float(elev) < el_max:
                        xs.append(x_idx)
                        ys.append(y_idx)
            if xs:
                pltx.scatter(xs, ys, color=color, marker='▪', label=label)

        pltx.yticks(list(range(n_ants)), antenna_names)
        pltx.ylim(-0.5, n_ants - 0.5)

        # Show three x-axis ticks: start, middle, end
        mid = n_times // 2
        pltx.xticks([0, mid, n_times - 1], [time_labels[0], time_labels[mid], time_labels[-1]])
        pltx.xlim(-legend_offset, n_times - 0.5)

        pltx.title(src_name)
        pltx.show()

    def _summary_tui(self):
        """Prints the observation summary to the terminal: band/time/duration, setup
        (data rate, bandwidth, subbands, channels, polarizations), the participating
        (and excluded) stations, the source list per scan block, and either a
        phase-referencing feasibility warning or (if no sources are defined) the
        expected thermal noise for a source at +/-45 degrees elevation.
        """
        rprint("\n[bold green]VLBI observation[/bold green]")
        rprint(f"To be conducted at {self.band.replace('cm', ' cm')} ", end='')
        if self.fixed_time:
            rprint(f"from {self.times[0].strftime('%d %b %Y %H:%M')}–"
                   f"{self.times[-1].strftime('%H:%M')} UTC")
        else:
            rprint("at unspecified times.")

        if self.duration is not None:
            rprint(f"With a total duration of {optimal_units(self.duration, [u.h, u.min, u.s]):.01f}.")

        rprint("\n[bold green]Setup[/bold green]")
        if None not in (self.datarate, self.bandwidth, self.subbands):
            val = optimal_units(self.datarate, [u.Gbit/u.s, u.Mbit/u.s])
            rprint(f"\nData rate of {val.value:.0f} {val.unit.to_string('unicode')}, "
                   f"producing a total bandwidth of {optimal_units(self.bandwidth, [u.MHz, u.GHz])}, "
                   f"divided in {self.subbands} x {int(self.bandwidth.value/self.subbands)}-"
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

        # Report stations that were requested but excluded
        excluded_msgs: list[str] = []
        for code, reason in self._excluded_stations.items():
            excluded_msgs.append(f"{code} ({reason})")
        if self.sources() and self.fixed_time:
            can_obs = self.can_be_observed()
            never_visible = set()
            for blk_obs in can_obs.values():
                for ant, visible in blk_obs.items():
                    if not visible:
                        never_visible.add(ant)
            # Only flag stations that cannot observe ANY block
            all_blocks_invisible = never_visible.copy()
            for blk_obs in can_obs.values():
                all_blocks_invisible &= {ant for ant, v in blk_obs.items() if not v}
            for ant in sorted(all_blocks_invisible):
                excluded_msgs.append(f"{ant} (source not visible)")
        if excluded_msgs:
            rprint(f"    [yellow]Not observing: {', '.join(excluded_msgs)}[/yellow]")

        rprint("\n[bold green]Sources[/bold green]:")
        if len(self.scans) > 0:
            for ablockname, ablock in self.scans.items():
                rprint(f"    - [dim]ScanBlock[/dim] '{ablockname}'")
                rprint('      ' +
                       '\n      '.join([s.name + ' [dim](' + s.coord.to_string('hmsdms') + ')[/dim]'
                                        for s in ablock.sources()]))
                if not ablock.sources(sources.SourceType.PHASECAL) and self.frequency > 80*u.GHz:
                    rprint("[red]Phase-referencing is not feasible anymore at this band.\n"
                           "Slewing times would be too short due to the coherent time.\n"
                           "A bright target is thus mandatory.[/red]")

        else:
            rprint("[dim]No sources defined.[/dim]")
            if self.fixed_time or self.duration is not None:
                rms = self.thermal_noise()
                if rms is not None:
                    rprint("\n[bold green]Expected outcome[/bold green]:")
                    val = optimal_units(rms, [u.Jy/u.beam, u.mJy/u.beam, u.uJy/u.beam])
                    rprint("[bold]Thermal rms noise (for a +/- 45° elevation source)[/bold]: ", end='')
                    rprint(f"{val.value:.01f} {val.unit.to_string('unicode')}")
                    rprint("[dim](for a +/- 45° elevation source)[/dim]")
                else:
                    rprint("\n[yellow]Cannot compute thermal noise (no stations with this band).[/yellow]")

        print('\n')

    def plot_visibility(self, gui: bool = True, tui: bool = True):
        """Shows visibility plots for the observation's sources.

        Unlike `summary()`, both variants are shown when both flags are True
        (they are not mutually exclusive here).

        Parameters
        ----------
        gui : bool, optional
            If True, shows the GUI visibility plot (via the `gui.plots` module). Default is True.
        tui : bool, optional
            If True, shows the terminal (TUI) visibility plot. Default is True.
        """
        if gui:
            self._plot_visibility_gui()

        if tui:
            self._plot_visibility_tui()

    def _plot_visibility_tui(self):
        """Show plots on the Terminal User Interface with the different sources and
        when they are visible within the observation.

        Displays per-source visibility bars (not per-block AND) so the user can see
        each source's observability independently.  The 'when everyone can observe'
        and 'optimal visibility range' messages refer to the TARGET source only.
        Sun proximity warnings identify the specific source that is too close.
        """
        if not self.fixed_time:
            rprint("[bold green]Searching for suitable GST range "
                   "(no pre-defined observing time)[/bold green]\n")
            self.times = Time('2025-09-21', scale='utc') + np.arange(0.0, 1.005, 0.01)*u.day
            doing_gst = True
        else:
            doing_gst = False

        per_src = self.per_source_observable()
        per_src_elev = self.per_source_elevations()
        rms_noise = self.thermal_noise()
        ontarget_time = self.ontarget_time
        sun_per_src = self.sun_constraint_per_source(times=self._REF_YEAR if doing_gst else None)
        sun_limit = self.sun_limiting_epochs()
        if doing_gst:
            gstimes = self.gstimes
            localtimes = self.times[:]
            self.times = None

        rprint("[bold]The blocks are observable for:[/bold]")
        for ablockname, ablock in self.scans.items():
            rprint(f"    - '{ablockname}':")
            # Show only science targets (TARGET or PULSAR) in the CLI plot; fall back to all if none exist.
            block_sources = (ablock.sources(sources.SourceType.TARGET) or
                             ablock.sources(sources.SourceType.PULSAR) or
                             ablock.sources())
            for src in block_sources:
                src_vis = per_src.get(src.name, {})
                src_elev = per_src_elev.get(src.name, {})
                if not src_vis:
                    continue
                # Check if source is always observable by all stations
                always_all = all(np.all(v) for v in src_vis.values())
                always_some = [ant for ant, v in src_vis.items() if np.all(v)]
                rprint(f"      [dim]{src.name}[/dim]", end='')
                if always_all:
                    rprint(" [dim](always observable)[/dim]")
                elif always_some:
                    cant = [ant for ant in self.stations.station_codenames if ant not in always_some]
                    rprint(f" [dim](always observable by everyone but {','.join(cant)})[/dim]")
                else:
                    rprint(' [dim](nobody can observe it all the time)[/dim]')

                # Build time-axis labels then hand off to plotext
                if doing_gst:
                    t_wrap = (24*u.hourangle
                              if np.abs(localtimes[-1].mjd - localtimes[0].mjd - 1) < 0.1
                              else 0.0*u.hourangle)
                    time_labels = [f"{gst.to_string(sep=':', fields=2, pad=True)} GST"
                                   for gst in gstimes]
                    time_labels[-1] = (gstimes[-1] + t_wrap).to_string(sep=':', fields=2, pad=True)
                else:
                    time_labels = [f"{t.strftime('%H:%M')} UTC" for t in self.times.datetime]

                self._plot_source_elevation_terminal(src.name, src_vis, src_elev, time_labels)

            # "When everyone can observe" and "optimal visibility" — use block-level
            if doing_gst:
                when_everyone = self.when_is_observable(mandatory_stations='all',
                                                        return_gst=True)[ablockname]
                if len(when_everyone) > 0:
                    rprint("\n[bold]All antennas can observe the block simultaneously at: [/bold]", end='')
                    rprint(', '.join([t1.to_string(sep=':', fields=2, pad=True) + '-' +
                           t2.to_string(sep=':', fields=2, pad=True) +
                           ' GST' for t1, t2 in when_everyone]))
                else:
                    rprint("\nThe block cannot be observed by all stations at the same time.")
            else:
                when_everyone = self.when_is_observable(mandatory_stations='all')[ablockname]
                if len(when_everyone) > 0:
                    rprint("\n[bold]All antennas can observe the block simultaneously at: [/bold]", end='')
                    rprint(', '.join([t1.strftime('%d %b %Y %H:%M')+'-'+t2.strftime('%H:%M') +
                           ' UTC' for t1, t2 in when_everyone]))
                else:
                    rprint("\nThe block cannot be observed by all stations at the same time.")

            min_stat = 3 if len(self.stations) > 3 else min(2, len(self.stations))
            if len(self.stations) > 2:
                rprint(f"[bold]Optimal visibility window (> {min_stat} antennas simultaneously):[/bold] ", end='')
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

            # Sun constraint — per source, so the user knows which source is the problem
            if doing_gst:
                if len(sun_limit[ablockname]) > 0:
                    offending = [(s.name, sun_per_src[s.name])
                                 for s in block_sources if sun_per_src.get(s.name) is not None]
                    src_label = ', '.join(s for s, _ in offending) if offending else 'a source in this block'
                    min_sep = min((sep for _, sep in offending), default=None) if offending else None
                    sep_str = f" (min separation {min_sep:.1f})" if min_sep is not None else ''
                    rprint(f"[bold red]Note that the Sun is too close to {src_label}{sep_str}[/bold red]",
                           end='')
                    t0, t1 = sun_limit[ablockname][0].datetime, sun_limit[ablockname][-1].datetime
                    if t0 == t1:
                        rprint(f"[bold red] on {t0.strftime('%d %b')}![/bold red]")
                    elif t0.month == t1.month:
                        rprint(f"[bold red] on {t0.day}-{t1.day} {t0.strftime('%b')}![/bold red]")
                    elif t0.year == t1.year:
                        rprint(f"[bold red] on {t0.strftime('%d %b')} to "
                               f"{t1.strftime('%d %b')}![/bold red]")
                    else:
                        rprint(f"[bold red]{t0.strftime('%d %b %Y')} to "
                               f"{t1.strftime('%d %b %Y')}![/bold red]")
            else:
                offending = [(s.name, sun_per_src[s.name])
                             for s in block_sources if sun_per_src.get(s.name) is not None]
                for src_name, sep in offending:
                    rprint(f"[bold red]Note that the Sun is too close to {src_name} during "
                           f"this observation (separation of {sep:.1f}).[/bold red]")

            rprint("[bold]Expected rms thermal noise for the target source: [/bold]", end='')
            if rms_noise is None:
                rprint("[yellow]Cannot be computed (not enough simultaneous antenna coverage).[/yellow]")
            else:
                for src, rms in rms_noise.items():
                    if rms is None:
                        if src in self.sourcenames_in_block(ablockname):
                            rprint(f"[yellow]{src}: cannot be computed "
                                   "(not enough simultaneous antenna coverage).[/yellow]")
                        continue
                    if any([s.type is sources.SourceType.TARGET for s in self.sources()]):
                        if src in self.sourcenames_in_block(ablockname, sources.SourceType.TARGET):
                            val = optimal_units(rms, [u.Jy/u.beam, u.mJy/u.beam, u.uJy/u.beam])
                            rprint(f"{src}: {val.value:.02f} {val.unit.to_string('unicode')}")
                            val = optimal_units(ontarget_time[src], [u.h, u.min, u.s, u.ms])
                            rprint("[dim]for a total on-source time of ~ "
                                   f"{val.value:.2f} {val.unit.to_string('unicode')} "
                                   f"(assuming the total observing time).[/dim]")
                    else:
                        if src in self.sourcenames_in_block(ablockname):
                            val = optimal_units(rms, [u.Jy/u.beam, u.mJy/u.beam, u.uJy/u.beam])
                            rprint(f"{src}: {val.value:.2f} {val.unit.to_string('unicode')}.")

            print('\n')

    def _plot_visibility_gui(self):
        """Show plots with the different sources and when they are visible within the
        observation
        """
        # Lazy import: plots pulls in plotly/matplotlib, only needed for GUI output.
        from vlbiplanobs.gui import plots

        if self.scans is None:
            # rprint("No scans have been defined.")
            sys.exit(0)

        figs = plots.elevation_plot(self)
        figs.show()

    def print_baseline_sensitivities(self):
        """Prints the sensitivity of all baselines in the observation per one minute of time
        """
        rprint("[bold green]Sensitivity per baseline (in mJy/beam)[/bold green]\n")
        rprint(" "*8, end='')
        for s in self.stations:
            rprint(f"{s.codename:6}", end='')

        rprint('\n' + '------' * (len(self.stations) + 1))
        for i in range(len(self.stations)):
            rprint(f"{self.stations[i].codename:3} | {' '*i*6}", end='')
            for j in range(i, len(self.stations)):
                try:
                    temp = self.baseline_sensitivity(self.stations[i].codename,
                                                     self.stations[j].codename).to(u.mJy/u.beam).value
                    rprint(f"{temp:6.3}", end='')
                except TypeError:
                    rprint("     ", end='')

            rprint('')
