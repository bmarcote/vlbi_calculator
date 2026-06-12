from typing import Optional
import base64
import json
import importlib.metadata
from packaging.version import Version
from urllib.parse import quote, unquote
from dash import html, Output, Input, State, callback, ctx, no_update, clientside_callback, ALL, MATCH
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
from astropy import coordinates as coord
from astropy import units as u
from furl import furl
from vlbiplanobs import freqsetups as fs
from vlbiplanobs import sources
from vlbiplanobs import observation
from vlbiplanobs import cli
from vlbiplanobs.gui import inputs, plots


@callback([Output('band-slider', 'marks'),
           Output({'type': 'network-label-band', 'index': ALL}, 'children'),
           Output({'type': 'badge-band-ant', 'index': ALL}, 'children')],
          Input('switch-band-label', 'value'))
def change_band_labels(show_wavelengths: bool):
    return {i: label for i, label in enumerate(inputs.pick_band_labels(show_wavelengths))}, \
           [inputs.network_band_labels(network, show_wavelengths)
            for network in observation._NETWORKS], \
           [inputs.print_table_bands_sefds(ant, show_wavelengths)
            for ant in observation._STATIONS]


@callback([Output({'type': 'network-switch', 'index': ALL}, 'disabled'),
           Output({'type': 'network-card', 'index': ALL}, 'style')],
          [Input('band-slider', 'value'),
           Input('store-targets', 'data')],
          State({'type': 'network-card', 'index': ALL}, 'style'))
def enable_networks_with_band(band_index: int, target_specs: Optional[list[str]],
                              card_styles):
    """Disable/enable network cards based on selected band and source observability.

    A network is enabled only when it supports the selected band (if any) AND any of the
    specified target sources is observable by at least 3 of its stations (if any sources
    are given).
    """
    the_band = inputs.band_from_index(band_index) if band_index != 0 else None
    band_ok = {k: (the_band in net.observing_bands if the_band else True)
               for k, net in observation._NETWORKS.items()}

    source_ok = {k: True for k in observation._NETWORKS}
    target_specs = target_specs or []
    if target_specs:
        # Use the first target spec only as a quick gating heuristic: full validation
        # happens on each target individually inside the compute callback.
        first_spec = next((s for s in target_specs if s and s.strip()), None)
        if first_spec:
            try:
                src = sources.Source.source_from_str(first_spec)
                for net_key, network in observation._NETWORKS.items():
                    if not band_ok[net_key]:
                        continue
                    check_band = the_band if the_band else list(network.observing_bands)[0]
                    try:
                        obs_obj = observation.Observation(
                            band=check_band, stations=network, times=None, duration=1 * u.hour,
                            datarate=network.max_datarate(check_band),
                            scans={src.name: sources.ScanBlock([sources.Scan(src, duration=5 * u.min)])}
                        )
                        observable = obs_obj.when_is_observable(min_stations=3)
                        source_ok[net_key] = len(observable[src.name]) > 0
                    except Exception:
                        source_ok[net_key] = False
            except Exception:
                pass

    enabled = {k: band_ok[k] and source_ok[k] for k in observation._NETWORKS}
    new_card_styles = tuple({k: v if k != 'opacity' else (1.0 if enabled[nk] else 0.2)
                             for k, v in style.items()}
                            for nk, style in zip(observation._NETWORKS, card_styles))
    return tuple(not enabled[k] for k in observation._NETWORKS), new_card_styles


@callback([Output('datarate', 'options'),
           Output('datarate', 'value', allow_duplicate=True),
           Output('channels', 'value', allow_duplicate=True),
           Output('subbands', 'value', allow_duplicate=True)],
          [Input('switch-specify-continuum', 'value'),
           Input('band-slider', 'value'),
           Input({'type': 'network-switch', 'index': ALL}, 'value')],
          [State('datarate', 'value'),
           State('store-prev-channels', 'data'),
           State('store-prev-subbands', 'data')],
          prevent_initial_call=True)
def prioritize_spectral_line(do_spectral_line: bool, band: int, network_bools: list[bool],
                             datarate: int = 2048, prev_channels: int = 64, prev_subbands: int = 8):
    if band == 0:
        raise PreventUpdate

    if not [None for nb, nn in zip(network_bools, observation._NETWORKS.values())
            if nb and nn.has_band(inputs.band_from_index(band))]:
        return no_update, no_update, no_update, no_update

    network_names = [nn for nb, nn in zip(network_bools, observation._NETWORKS) if nb]
    try:
        the_band = inputs.band_from_index(band)
        if 'EVN' in network_names and observation._NETWORKS['EVN'].has_band(the_band):
            dr = observation._NETWORKS['EVN'].max_datarate(the_band)
            max_datarate = int(dr.value) if dr is not None else 2048
        else:
            rates = [observation._NETWORKS[net].max_datarate(the_band)
                     for net in network_names if observation._NETWORKS[net].has_band(the_band)]
            valid_rates = [r.value for r in rates if r is not None]
            max_datarate = int(min(valid_rates)) if valid_rates else 2048
    except (AttributeError, ValueError):
        raise PreventUpdate

    # This is here to check if this avoids the state when somehow datarate value is None, which I do not see why
    if datarate is None:
        datarate = 2048
    return  \
        tuple({'value': dr, 'label': html.Span([drl], style={'color': '#888888'
                                                             if dr > max_datarate else '#000000'})}
              for dr, drl in fs.data_rates.items()), \
        32 if do_spectral_line else max_datarate if max_datarate is not None else datarate, \
        4096 if do_spectral_line else prev_channels, \
        1 if do_spectral_line else prev_subbands


@callback(Output({'type': 'antenna-chip', 'index': ALL}, 'disabled'),
          [Input('band-slider', 'value'),
           Input('switch-specify-e-evn', 'value')])
def enable_antennas_with_band(band_index: int, do_e_evn: bool):
    grouped_codenames = {s.codename for slist in inputs.station_groups().values() for s in slist}
    ungrouped = [ant for ant in observation._STATIONS if ant.codename not in grouped_codenames]
    if band_index == 0:
        return [False for _ in ungrouped]

    return [not ant.has_band(inputs.band_from_index(band_index)) or (do_e_evn and not ant.real_time)
            for ant in ungrouped]


@callback(
    [Output({'type': 'group-is-selected', 'index': MATCH}, 'data'),
     Output('switches-antennas', 'value', allow_duplicate=True)],
    Input({'type': 'group-toggle-btn', 'index': MATCH}, 'n_clicks'),
    [State({'type': 'group-is-selected', 'index': MATCH}, 'data'),
     State({'type': 'group-active-codename', 'index': MATCH}, 'data'),
     State('switches-antennas', 'value')],
    prevent_initial_call=True
)
def toggle_group_chip(n_clicks, is_selected, active_codename, current_antennas):
    """Toggle a grouped antenna chip on/off, adding or removing its active codename
    from the switches-antennas selection."""
    if n_clicks is None:
        raise PreventUpdate
    new_selected = not bool(is_selected)
    current_set = set(current_antennas or [])
    if new_selected:
        current_set.add(active_codename)
    else:
        current_set.discard(active_codename)
    return new_selected, list(current_set)


@callback(
    [Output({'type': 'group-active-codename', 'index': ALL}, 'data', allow_duplicate=True),
     Output({'type': 'group-is-selected', 'index': ALL}, 'data', allow_duplicate=True),
     Output('switches-antennas', 'value', allow_duplicate=True)],
    Input({'type': 'group-menu-item', 'index': ALL}, 'n_clicks'),
    [State({'type': 'group-active-codename', 'index': ALL}, 'data'),
     State({'type': 'group-is-selected', 'index': ALL}, 'data'),
     State('switches-antennas', 'value')],
    prevent_initial_call=True
)
def switch_group_config(menu_clicks, active_codenames, is_selected_list, current_antennas):
    """Switch the active configuration for a grouped antenna when the user picks a menu item.
    Always selects the group (turns it on) when a config is chosen — picking a config implies
    including that antenna in the observation. Removes the old codename and adds the new one
    to switches-antennas.
    """
    if not ctx.triggered_id:
        raise PreventUpdate
    triggered_index = ctx.triggered_id.get('index', '')
    if '__' not in triggered_index:
        raise PreventUpdate
    group_name_triggered, new_codename = triggered_index.split('__', 1)

    # Build ordered list of group names matching the ALL store order
    group_names = [item['id']['index'] for item in ctx.states_list[0]]
    try:
        group_idx = group_names.index(group_name_triggered)
    except ValueError:
        raise PreventUpdate

    current_codename = active_codenames[group_idx]

    new_active_codenames = list(active_codenames)
    new_active_codenames[group_idx] = new_codename

    new_is_selected = [no_update] * len(is_selected_list)
    new_is_selected[group_idx] = True

    current_set = set(current_antennas or [])
    current_set.discard(current_codename)
    current_set.add(new_codename)

    return new_active_codenames, new_is_selected, list(current_set)


@callback(
    Output({'type': 'group-menu-item', 'index': ALL}, 'className'),
    [Input({'type': 'group-active-codename', 'index': ALL}, 'data'),
     Input('switches-antennas', 'value')],
)
def highlight_active_menu_item(active_codenames, _switches):
    """Set className on each menu item to mark the active configuration.
    Uses a CSS class rather than Mantine style props so the highlight survives
    Dash re-renders triggered by unrelated callbacks (e.g. observation compute).
    Re-triggered by switches-antennas changes to re-apply after any full update.
    """
    active_set = set(active_codenames)

    classes = []
    for item in ctx.outputs_list[0]:
        index = item['id']['index']  # format: 'groupname__codename'
        codename = index.split('__', 1)[1] if '__' in index else ''
        classes.append('group-menu-item-active' if codename in active_set else '')
    return classes


@callback(
    [Output({'type': 'group-toggle-btn', 'index': MATCH}, 'children'),
     Output({'type': 'group-toggle-btn', 'index': MATCH}, 'className'),
     Output({'type': 'group-chip-wrapper', 'index': MATCH}, 'style')],
    [Input({'type': 'group-active-codename', 'index': MATCH}, 'data'),
     Input({'type': 'group-is-selected', 'index': MATCH}, 'data'),
     Input('band-slider', 'value'),
     Input('switch-specify-e-evn', 'value')],
)
def update_group_chip_appearance(active_codename, is_selected, band_index, do_e_evn):
    """Update the toggle button label and styling to reflect the active config and
    selected/disabled state. Also dims the group chip wrapper when the active config
    cannot observe the selected band."""
    ant = observation._STATIONS[active_codename]
    group_name = ctx.outputs_list[0]['id']['index']
    label = group_name.upper()

    disabled = False
    if band_index and band_index != 0:
        band = inputs.band_from_index(band_index)
        disabled = not ant.has_band(band) or (do_e_evn and not ant.real_time)

    css_class = 'btn-group-chip-toggle '
    if disabled:
        css_class += 'btn-group-chip-disabled'
    elif is_selected:
        css_class += 'btn-group-chip-on'
    else:
        css_class += 'btn-group-chip-off'

    wrapper_style = {'display': 'inline-flex', 'align-items': 'center', 'align-self': 'center',
                     'opacity': '0.45' if disabled else '1.0'}
    return label, css_class, wrapper_style


@callback(
    [Output('switches-antennas', 'value', allow_duplicate=True)] +
    [Output({'type': 'group-is-selected', 'index': gname}, 'data', allow_duplicate=True)
     for gname in inputs.station_groups()],
    Input({'type': 'network-switch', 'index': ALL}, 'value'),
    [State('switches-antennas', 'value')] +
    [State({'type': 'group-active-codename', 'index': gname}, 'data')
     for gname in inputs.station_groups()],
    prevent_initial_call=True
)
def update_selected_antennas_from_networks(networks, current_antennas, *group_active_codenames):
    """When a network switch is toggled, update both the ungrouped chip selection and
    the is-selected state for each grouped antenna whose active codename is in the network."""
    current_antennas = set(current_antennas)
    ants2include = set()
    ants2exclude = set()

    for selected, network in zip(networks, observation._NETWORKS.values()):
        if selected:
            ants2include.update(network.station_codenames)
        else:
            ants2exclude.update(network.station_codenames)

    ants2exclude -= ants2include
    new_antennas = current_antennas - ants2exclude | ants2include

    group_selected_states = []
    for active_codename in group_active_codenames:
        if active_codename in ants2include:
            group_selected_states.append(True)
        elif active_codename in ants2exclude:
            group_selected_states.append(False)
        else:
            # No network touched this group's active config — keep current state
            group_selected_states.append(no_update)

    return [list(new_antennas)] + group_selected_states


# Epoch toggle: show/hide the date+time fields based on the switch.
# This MUST be a serverside callback (not clientside): on page load the serverside
# initial call receives the persistence-restored switch value, so a previously-ticked
# switch correctly reveals the epoch fields. A clientside callback races against
# persistence restoration and fires with the stale default (False), leaving the fields
# hidden until the user toggles the switch off and on again.
@callback(Output('epoch-selection-div', 'className'),
          Input('switch-specify-epoch', 'value'))
def toggle_epoch_selection(value: bool) -> str:
    """Return '' to show the epoch date/time fields when the switch is on, else 'd-none'."""
    return '' if value else 'd-none'


# --------------------------------------------------------------------------------------
# Multi-target source management
# --------------------------------------------------------------------------------------
# The list of target source specs (strings) lives in `dcc.Store(id='store-targets')`.
# The user manages it through the modal opened from the source-and-epoch panel.

clientside_callback(
    "function(n_open, n_close, is_open) { return !is_open; }",
    Output('source-modal', 'is_open'),
    [Input('button-open-source-modal', 'n_clicks'),
     Input('button-close-source-modal', 'n_clicks')],
    State('source-modal', 'is_open'),
    prevent_initial_call=True
)


def _validate_source_spec(spec: str) -> tuple[bool, str]:
    """Validate a single source spec string.

    Returns (ok, message). The message describes the resolved coordinates on success
    or the reason for failure otherwise.
    """
    spec = (spec or '').strip()
    if not spec:
        return False, "Empty source name/coordinates."
    if len(spec) > 80:
        return False, "Name too long."
    try:
        src = sources.Source.source_from_str(spec)
    except ValueError:
        return False, "Wrong coordinates format. Use 'hh:mm:ss dd:mm:ss' or 'XXhXXmXXs XXdXXmXXs'."
    except coord.name_resolve.NameResolveError:
        return False, "Unrecognized source name (not found in SIMBAD/NED/VizieR)."
    except Exception as exc:
        return False, f"Could not parse source: {exc}"
    return True, src.coord.to_string('hmsdms', precision=3)


@callback([Output('store-targets', 'data'),
           Output('modal-source-input', 'value'),
           Output('modal-source-feedback', 'children'),
           Output('modal-source-feedback', 'className'),
           Output('upload-sources-feedback', 'children'),
           Output('upload-sources-feedback', 'className')],
          [Input('button-add-source', 'n_clicks'),
           Input('modal-source-input', 'n_submit'),
           Input('upload-sources', 'contents'),
           Input('button-clear-sources', 'n_clicks'),
           Input({'type': 'modal-remove-target', 'index': ALL}, 'n_clicks')],
          [State('modal-source-input', 'value'),
           State('upload-sources', 'filename'),
           State('store-targets', 'data')],
          prevent_initial_call=True)
def manage_target_sources(add_clicks, submit_n, upload_contents, clear_clicks,
                          remove_clicks, modal_input, upload_filename, store_data):
    """Single dispatcher for all target-source mutations driven from the modal.

    Handles: add via text input, add via file upload, remove individual, and clear all.
    Returns the updated store plus user feedback messages.
    """
    targets: list[str] = list(store_data or [])
    triggered = ctx.triggered_id
    feedback_text, feedback_class = no_update, no_update
    upload_text, upload_class = no_update, no_update
    new_modal_input = no_update

    def _append_unique(spec: str) -> bool:
        spec = spec.strip()
        if spec and spec not in targets:
            targets.append(spec)
            return True
        return False

    if triggered == 'button-clear-sources':
        targets = []
        feedback_text, feedback_class = "All sources cleared.", 'form-text text-warning'

    elif triggered in ('button-add-source', 'modal-source-input'):
        ok, msg = _validate_source_spec(modal_input or '')
        if not ok:
            feedback_text, feedback_class = msg, 'form-text text-danger'
        elif (modal_input or '').strip() in targets:
            feedback_text = f"'{modal_input.strip()}' is already in the list."
            feedback_class = 'form-text text-warning'
        else:
            _append_unique(modal_input or '')
            new_modal_input = ''
            feedback_text = f"Added (resolved to {msg})."
            feedback_class = 'form-text text-success'

    elif triggered == 'upload-sources':
        if upload_contents is None:
            raise PreventUpdate
        try:
            _, content_string = upload_contents.split(',', 1)
            decoded = base64.b64decode(content_string).decode('utf-8', errors='replace')
        except Exception as exc:
            upload_text = f"Could not read file: {exc}"
            upload_class = 'form-text text-danger'
        else:
            added, skipped, invalid = 0, 0, []
            for raw in decoded.splitlines():
                line = raw.strip()
                if not line or line.startswith('#'):
                    continue
                ok, msg = _validate_source_spec(line)
                if not ok:
                    invalid.append(f"{line!r}: {msg}")
                    continue
                if _append_unique(line):
                    added += 1
                else:
                    skipped += 1
            label = upload_filename or 'file'
            parts = [f"Loaded {added} source(s) from {label}."]
            if skipped:
                parts.append(f"{skipped} already present.")
            if invalid:
                parts.append(f"{len(invalid)} invalid line(s) ignored.")
            upload_text = ' '.join(parts)
            upload_class = 'form-text text-success' if added else 'form-text text-warning'

    elif isinstance(triggered, dict) and triggered.get('type') == 'modal-remove-target':
        idx = triggered.get('index')
        # n_clicks may fire on initial render; ignore None entries.
        if 0 <= idx < len(targets) and any(remove_clicks or []):
            removed = targets.pop(idx)
            feedback_text = f"Removed '{removed}'."
            feedback_class = 'form-text text-muted'
        else:
            raise PreventUpdate
    else:
        raise PreventUpdate

    return targets, new_modal_input, feedback_text, feedback_class, upload_text, upload_class


@callback([Output('target-chips-display', 'children'),
           Output('modal-source-list', 'children')],
          Input('store-targets', 'data'))
def render_target_sources(targets: Optional[list[str]]):
    """Render the chip list shown in the inputs panel and the detailed list inside the modal.

    Each modal entry has an X button (pattern-matching id `modal-remove-target`).
    """
    targets = targets or []

    if not targets:
        chips = html.Small("No target sources added yet.", className='text-muted')
        modal_list = html.P("No sources added yet.",
                            className='text-muted text-center my-3')
        return chips, modal_list

    chip_children = []
    for spec in targets:
        chip_children.append(
            dbc.Badge(spec, pill=True, className='me-1 mb-1',
                      style={'font-size': '0.85em', 'padding': '0.4em 0.7em',
                             'background-color': '#004990', 'color': 'white'}))
    chips = html.Div(chip_children, style={'max-height': '6em', 'overflow-y': 'auto'})

    modal_items = []
    for i, spec in enumerate(targets):
        modal_items.append(
            html.Li(className='list-group-item d-flex justify-content-between align-items-center',
                    children=[
                        html.Span(spec, className='text-break'),
                        dbc.Button(html.I(className='fa fa-times'),
                                   id={'type': 'modal-remove-target', 'index': i},
                                   color='danger', outline=True, size='sm',
                                   n_clicks=0,
                                   title=f"Remove {spec}")
                    ]))
    modal_list = html.Ul(modal_items, className='list-group list-group-flush')

    return chips, modal_list


@callback(Output('bandwidth-label', 'children'),
          [Input('datarate', 'value'),
           Input('pols', 'value'),
           Input('channels', 'value'),
           Input('subbands', 'value')],
          State('switch-specify-continuum', 'value'))
def update_bandwidth_label(datarate: int, npols: int, chans: int, subbands: int, do_spectral_line: bool):
    """Updates the total bandwidth label as a function of the selected datarate and number of
    polarizations. Returns a string with the value and units.
    """
    if None in (datarate, npols, chans, subbands):
        raise PreventUpdate

    # Either 1 or 2 pols per station:
    return "The maximum bandwidth is " \
           f"{cli.optimal_units(int(datarate)*u.MHz/((npols % 3 + npols//3)*2*2), [u.GHz, u.MHz, u.kHz]):.0f}."


clientside_callback(
    "function(n_clicks, is_open) { return n_clicks ? !is_open : dash_clientside.no_update; }",
    Output("more-info-modal", "is_open"),
    Input("more-info-button", "n_clicks"),
    State("more-info-modal", "is_open"),
    prevent_initial_call=True
)


# Per-tab UV antenna-highlight callback. Each tab embeds:
#   - a `dcc.Dropdown(id={'type':'select-ant-uv','index':<src>})`
#   - a `dcc.Graph(id={'type':'fig-uv','index':<src>})`
#   - a `dcc.Store(id={'type':'store-uv-data','index':<src>}, data=<serialized uv>)`
# This pattern-matching callback updates the figure when the user picks antennas.
@callback(Output({'type': 'fig-uv', 'index': MATCH}, 'figure'),
          Input({'type': 'select-ant-uv', 'index': MATCH}, 'value'),
          State({'type': 'store-uv-data', 'index': MATCH}, 'data'),
          prevent_initial_call=True)
def update_uv_figure(highlight_antennas: list[str], uv_data: dict):
    if uv_data is None:
        raise PreventUpdate
    return plots.uvplot_from_data(uv_data, highlight_antennas)


@callback([Output('error_duration', 'children'),
           Output('error_duration', 'className'),
           Output('duration', 'className')],
          Input('duration', 'value'), prevent_initial_call=True)
def check_initial_obstime(duration: Optional[int | float]):
    """Verify the introduced times/dates for correct values.
    Once the user has introduced all values for the start and end of the observation,
    it guarantees that they have the correct shape:
        - the duration of the observation is > 0 hours.
        - The total observing length is less than five days (value chosen for computational reasons).
    """
    if duration is None:
        return no_update, no_update, no_update

    if (not isinstance(duration, float)) and (not isinstance(duration, int)):
        return 'Must be a number',  'form-text text-danger', 'form-control is-invalid'

    if duration < 0.1:
        return 'Duration too short. Minimum duration is 0.1 h, but you can check the instantaneous ' \
               'sensitivity in the provided values (they are also calculated per minute integration).', \
               'form-text text-danger', 'form-control is-invalid'
    elif duration > 4*24:
        return 'Must be shorter than 4 days',  'form-text text-danger', 'form-control is-invalid'

    return "", no_update, 'form-control'


# Theme toggle callback
clientside_callback(
    """
    function(n_clicks) {
        if (!n_clicks) {
            return dash_clientside.no_update;
        }
        if (window.planobsTheme) {
            const newTheme = window.planobsTheme.toggle();
            return newTheme === 'dark' ? 'Switch to Light Mode' : 'Switch to Dark Mode';
        }
        return dash_clientside.no_update;
    }
    """,
    Output("theme-toggle-text", "children"),
    Input("theme-toggle-btn", "n_clicks"),
    prevent_initial_call=True
)

# list of component IDs that constitute the user's configuration
export_component_ids = [
    'switch-band-label',
    'band-slider',
    'duration',
    'onsourcetime',
] + [
    {'type': 'network-switch', 'index': network_name}
    for network_name in observation._NETWORKS
] + [
    'switches-antennas',
    'switch-specify-epoch',
    'startdate',
    'starttime',
    'switch-specify-e-evn',
    'switch-specify-continuum',
    'datarate',
    'subbands',
    'channels',
    'pols',
    'inttime',
]

current_version = Version(importlib.metadata.version('vlbiplanobs'))

json_config = '[' + ','.join(f'[{json.dumps(export_component_ids[i])}, args[{i}]]' for i in range(len(export_component_ids))) + ']'
callback_javascript = f"""
    function(n_clicks, ...args) {{
        const value = '?targetversion={quote(str(current_version))}&config=' + encodeURIComponent(JSON.stringify({json_config}));
        if (window.opener && !window.opener.closed) {{
            window.opener.postMessage(value, '*'); // FIX set the true targetOrigin
            return [true, "Configuration sent to Polaris"];
        }}
        else {{
            navigator.clipboard.writeText(value);
            return [true, "Configuration copied to clipboard, paste in Polaris"];
        }}
    }}
"""
clientside_callback(
    callback_javascript,
    Output('export-alert', 'is_open'),
    Output('export-alert', 'children'),
    Input('export-state-of-the-system', 'n_clicks'),
    [State(id_, 'value') for id_ in export_component_ids],
    prevent_initial_call=True
    )

@callback(
    [Output(id_, 'value') for id_ in export_component_ids],
    # delete targetversion and config from the url parameters after parsing it
    Output('url', 'href'),
    Input('url', 'href')
    )
def url_open(href):
    parsed_href = furl(href)
    target_version = parsed_href.args.get('targetversion')
    config = parsed_href.args.get('config')
    if config is None:
        raise PreventUpdate
    if target_version is not None:
        target_version = Version(unquote(target_version))
        if target_version > current_version:
            # current running version too old??
            raise PreventUpdate
    config_list = json.loads(unquote(config))
    id_list, value_list = map(list, zip(*config_list))
    update_list = []
    for component_id in export_component_ids:
        try:
            index = id_list.index(component_id)
            update_list.append(value_list[index])
        except ValueError:
            update_list.append(no_update)
    if target_version is not None:
        del parsed_href.args['targetversion']
    if config is not None:
        del parsed_href.args['config']
    return update_list + [parsed_href.url]
