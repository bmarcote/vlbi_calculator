from typing import Optional
import json
import importlib.metadata
from packaging.version import Version
from urllib.parse import quote, unquote
from dash import html, Output, Input, State, callback, no_update, clientside_callback, ALL
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
          Input('band-slider', 'value'),
          State({'type': 'network-card', 'index': ALL}, 'style'))
def enable_networks_with_band(band_index: int, card_styles):
    if band_index == 0:
        return tuple(False for _ in observation._NETWORKS), \
                     tuple({k: v if k != 'opacity' else 1.0 for k, v in card_style.items()}
                           for card_style in card_styles)

    def opacity(x: bool) -> float:
        return 1.0 if x else 0.2

    new_card_styles: tuple = tuple({k: v if k != 'opacity' else
                                   opacity(inputs.band_from_index(band_index) in network.observing_bands)
                                   for k, v in card_style.items()}
                                   for network, card_style in zip(observation._NETWORKS.values(),
                                                                  card_styles))

    return tuple(inputs.band_from_index(band_index) not in network.observing_bands
                 for network in observation._NETWORKS.values()), new_card_styles


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
    if band_index == 0:
        return [False for _ in observation._STATIONS]

    return [not ant.has_band(inputs.band_from_index(band_index)) or (do_e_evn and not ant.real_time)
            for ant in observation._STATIONS]


@callback(Output('switches-antennas', 'value', allow_duplicate=True),
          Input({'type': 'network-switch', 'index': ALL}, 'value'),
          State('switches-antennas', 'value'),
          prevent_initial_call=True)
def update_selected_antennas_from_networks(networks, current_antennas):
    current_antennas = set(current_antennas)
    ants2include = set()
    ants2exclude = set()

    for selected, network in zip(networks, observation._NETWORKS.values()):
        if selected:
            ants2include.update(network.station_codenames)
        else:
            ants2exclude.update(network.station_codenames)

    # Exclude antennas that are not in the selected networks
    ants2exclude -= ants2include
    return tuple(current_antennas - ants2exclude | ants2include)


# Clientside callbacks for faster UI toggle responses
clientside_callback(
    "function(value) { return !value; }",
    Output('source-selection-div', 'hidden'),
    Input('switch-specify-source', 'value')
)

clientside_callback(
    "function(value) { return !value; }",
    Output('epoch-selection-div', 'hidden'),
    Input('switch-specify-epoch', 'value')
)


@callback([Output('error_source', 'children'),
           Output('error_source', 'className')],
          Input('source-input', 'value'))
def get_initial_source(source_coord: str):
    """Verifies that the introduced source coordinates have a right format.
    If they are correct, it does nothing. If they are incorrect, it shows an error label.
    """
    if (source_coord != 'hh:mm:ss dd:mm:ss') and (source_coord is not None) and (source_coord != ''):
        if len(source_coord) > 40:
            return "Name too long.", 'form-text text-danger'
        try:
            src = sources.Source.source_from_str(source_coord)
            return src.coord.to_string('hmsdms', precision=3), 'form-text text-success'
        except ValueError:
            return "Wrong coordinates.", 'form-text text-danger'
        except coord.name_resolve.NameResolveError:
            return "Unrecognized name. Use 'hh:mm:ss dd:mm:ss' or 'XXhXXmXXs XXdXXmXXs'", \
                   'form-text text-danger'

    return '', no_update


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
    Output("sensitivity-baseline-modal", "is_open"),
    Input("button-sensitivity-baseline", "n_clicks"),
    State("sensitivity-baseline-modal", "is_open"),
    prevent_initial_call=True
)

clientside_callback(
    "function(n_clicks, is_open) { return n_clicks ? !is_open : dash_clientside.no_update; }",
    Output("more-info-modal", "is_open"),
    Input("more-info-button", "n_clicks"),
    State("more-info-modal", "is_open"),
    prevent_initial_call=True
)


@callback([Output('fig-elevations', 'figure'),
           Output('fig-elevations2', 'figure')],
          Input('store-elev-data', 'data'), prevent_initial_call=True)
def render_elevation_plots(elev_data: dict):
    """Render elevation plots progressively from serialized data."""
    if elev_data is None:
        raise PreventUpdate
    return plots.elevation_plot_from_data(elev_data), plots.elevation_curves_from_data(elev_data)


@callback(Output('fig-uv-coverage', 'figure'),
          Input('store-uv-data', 'data'), prevent_initial_call=True)
def render_uv_plot(uv_data: dict):
    """Render UV coverage plot progressively from serialized data."""
    if uv_data is None:
        raise PreventUpdate
    return plots.uvplot_from_data(uv_data)


@callback(Output('fig-uv-coverage', 'figure', allow_duplicate=True),
          Input('select-antenna-uv-plot', 'value'),
          State('store-uv-data', 'data'), prevent_initial_call='initial_duplicate')
def update_uv_figure(highlight_antennas: list[str], uv_data: dict):
    if uv_data is None:
        raise PreventUpdate
    return plots.uvplot_from_data(uv_data, highlight_antennas)


@callback(Output('fig-worldmap', 'figure'),
          Input('store-worldmap-data', 'data'), prevent_initial_call=True)
def render_worldmap(worldmap_data: dict):
    """Render worldmap progressively from serialized data."""
    if worldmap_data is None:
        raise PreventUpdate
    return plots.worldmap_from_data(worldmap_data)


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
    'switch-specify-source',
    'source-input',
    'switch-specify-e-evn',
    'switch-specify-continuum',
    'datarate',
    'subbands',
    'channels',
    'pols',
    'inttime',
]

current_version = Version(importlib.metadata.version('vlbiplanobs'))

@callback(
    Output('url-store', 'data'),
    Input('export-state-of-the-system', 'n_clicks'),
    [State(id_, 'value') for id_ in export_component_ids],
    prevent_initial_call=True
    )
def clicked_export(n_clicks, *args):
    if not n_clicks:
        return no_update
    config = quote(json.dumps(
        [(id_, value) for id_, value in zip(export_component_ids, args)]))
    search = f'?targetversion={quote(str(current_version))}&config={config}'
    return search

clientside_callback(
    """
    function(value) {
        if (window.opener && !window.opener.closed) {
            window.opener.postMessage(value, '*'); // FIX set the true targetOrigin
            return [true, "Configuration sent to Polaris"];
        }
        else {
            navigator.clipboard.writeText(value);
            return [true, "Configuration copied to clipboard, paste in Polaris"];
        }
    }
    """,
    Output('export-alert', 'is_open'),
    Output('export-alert', 'children'),
    Input('url-store', 'data'),
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
    if (target_version is None) or (config is None):
        raise PreventUpdate
    target_version = Version(unquote(target_version))
    if target_version > current_version:
        # current running version too old??
        raise PreventUpdate
    config_list = json.loads(unquote(config))
    id_list, value_list = map(list, zip(*config_list))
    if id_list != export_component_ids:
        # mismatch in expected IDs, since only one version is supported,
        # this is currently an input error
        raise PreventUpdate
    del parsed_href.args['targetversion']
    del parsed_href.args['config']
    return value_list + [parsed_href.url]
