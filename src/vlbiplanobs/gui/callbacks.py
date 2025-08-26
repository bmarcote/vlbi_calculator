from typing import Optional
from dash import html, Output, Input, State, callback, no_update, Patch
from dash.exceptions import PreventUpdate
from astropy import coordinates as coord
from astropy import units as u
from vlbiplanobs import freqsetups as fs
from vlbiplanobs import sources
from vlbiplanobs import observation
from vlbiplanobs import cli
from vlbiplanobs.gui import inputs
from vlbiplanobs.gui import main as gui_main


@callback([Output('band-slider', 'marks'),
           [Output(f"network-{network}-label-band", 'children')
            for network in observation._NETWORKS],
           [Output(f"badge-band-ant-{ant.codename}", 'children')
            for ant in observation._STATIONS]],
          Input('switch-band-label', 'value'))
def change_band_labels(show_wavelengths: bool):
    return {i: label for i, label in enumerate(inputs.pick_band_labels(show_wavelengths))}, \
           tuple(inputs.network_band_labels(network, show_wavelengths)
                 for network in observation._NETWORKS), \
           tuple(inputs.print_table_bands_sefds(ant, show_wavelengths)
                 for ant in observation._STATIONS)


# @callback(Output("welcome-modal-shown", "data"),
#           Input("close-modal", "n_clicks"),
#           State("welcome-modal-shown", "data"))
# def update_modal_shown(n_clicks, modal_shown):
#     if n_clicks:
#         return True
#
#     return no_update
#
#
# @callback(Output("welcome-modal", "is_open"),
#           Input("close-modal", "n_clicks"),
#           [State("welcome-modal-shown", "data"),
#            State("welcome-modal", "is_open")])
# def toggle_modal(n_clicks, modal_shown, is_open):
#     if modal_shown is None:
#         return True
#
#     if n_clicks:
#         return not is_open
#
#     raise PreventUpdate


@callback([[Output(f"network-{network}", 'disabled') for network in observation._NETWORKS],
           [Output(f"network-{network}-card", 'style') for network in observation._NETWORKS]],
          Input('band-slider', 'value'),
          [State(f"network-{network}-card", 'style') for network in observation._NETWORKS])
def enable_networks_with_band(band_index: int, *card_styles):
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
           Output('datarate', 'value'),
           Output('channels', 'value'),
           Output('subbands', 'value')],
          [Input('switch-specify-continuum', 'value'),
           Input('band-slider', 'value'),
           [Input(f"network-{network}", 'value') for network in observation._NETWORKS]],
          State('datarate', 'value'))
def prioritize_spectral_line(do_spectral_line: bool, band: int, network_bools: list[bool], datarate: int = 2048):
    if band == 0:
        raise PreventUpdate

    if not [None for nb, nn in zip(network_bools, observation._NETWORKS.values())
            if nb and nn.has_band(inputs.band_from_index(band))]:
        return no_update, no_update, no_update, no_update

    network_names = [nn for nb, nn in zip(network_bools, observation._NETWORKS) if nb]
    try:
        # max_datarate is never None in the networks
        if 'EVN' in network_names and observation._NETWORKS['EVN'].has_band(inputs.band_from_index(band)):
            max_datarate = observation._NETWORKS['EVN'].max_datarate(inputs.band_from_index(band)).value
        else:
            max_datarate = int(min([observation._NETWORKS[net].max_datarate(inputs.band_from_index(band)).value
                                    for net in network_names
                                    if observation._NETWORKS[net].has_band(inputs.band_from_index(band))]))
    except AttributeError:
        raise PreventUpdate

    return  \
        tuple({'value': dr, 'label': html.Span([drl], style={'color': '#888888'
                                                             if dr > max_datarate else '#000000'})}
              for dr, drl in fs.data_rates.items()), \
        32 if do_spectral_line else max_datarate if max_datarate is not None else datarate, \
        4096 if do_spectral_line else gui_main._main_obs.prev_channels, \
        1 if do_spectral_line else gui_main._main_obs.prev_subbands


@callback([Output(f"chip-{ant.codename}", "disabled") for ant in observation._STATIONS],
          [Input('band-slider', 'value'),
           Input('switch-specify-e-evn', 'value')])
def enable_antennas_with_band(band_index: int, do_e_evn: bool):
    if band_index == 0:
        return tuple(False for _ in observation._STATIONS)

    return tuple(not ant.has_band(inputs.band_from_index(band_index)) or (do_e_evn and not ant.real_time)
                 for ant in observation._STATIONS)


@callback(Output('switches-antennas', 'value'),
          [Input(f"network-{network}", 'value') for network in observation._NETWORKS],
          State('switches-antennas', 'value'))
def update_selected_antennas_from_networks(*args):
    networks, current_antennas = args[:-1], set(args[-1])
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


@callback(Output('source-selection-div', 'hidden'),
          Input('switch-specify-source', 'value'))
def toggle_source_field(pick_source: bool):
    return not pick_source


@callback(Output('epoch-selection-div', 'hidden'),
          Input('switch-specify-epoch', 'value'))
def toggle_epoch_field(define_epoch: bool):
    return not define_epoch


@callback([Output('error_source', 'children'),
           Output('error_source', 'className'),
           Output('source-input', 'className')],
          Input('source-input', 'value'))
def get_initial_source(source_coord: str):
    """Verifies that the introduced source coordinates have a right format.
    If they are correct, it does nothing. If they are incorrect, it shows an error labinputs.
    """
    if (source_coord != 'hh:mm:ss dd:mm:ss') and (source_coord is not None) and (source_coord != ''):
        if len(source_coord) > 40:
            return "Name too long.", 'form-text text-danger', 'form-control'
        try:
            src = sources.Source.source_from_str(source_coord)
            return src.coord.to_string('hmsdms', precision=3), 'form-text text-success', 'form-control is-valid'
        except ValueError:
            return "Wrong coordinates.", 'form-text text-danger', 'form-control is-invalid'
        except coord.name_resolve.NameResolveError:
            return "Unrecognized name. Use 'hh:mm:ss dd:mm:ss' or 'XXhXXmXXs XXdXXmXXs'", \
                   'form-text text-danger', 'form-control is-invalid'

    return '', no_update, no_update


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

    if not do_spectral_line:
        gui_main._main_obs.prev_datarate = int(datarate)
        gui_main._main_obs.prev_channels = chans
        gui_main._main_obs.prev_subbands = subbands

    # Either 1 or 2 pols per station:
    return "The maximum bandwidth is " \
           f"{cli.optimal_units(int(datarate)*u.MHz/((npols % 3 + npols//3)*2*2), [u.GHz, u.MHz, u.kHz]):.0f}."


@callback(Output("sensitivity-baseline-modal", "is_open"),
          Input("button-sensitivity-baseline", "n_clicks"),
          State("sensitivity-baseline-modal", "is_open"), suppress_callback_exceptions=True)
def toggle_modal_baseline_sensitivity(n_clicks, is_open):
    if n_clicks is None:
        return no_update

    return not is_open


@callback(Output("more-info-modal", "is_open"),
          Input("more-info-button", "n_clicks"),
          State("more-info-modal", "is_open"), prevent_initial_call=True)
def toggle_modal_info_button(n_clicks, is_open):
    if n_clicks is None:
        return no_update

    return not is_open


@callback(Output('fig-uv-coverage', 'figure', allow_duplicate=True),
          Input('select-antenna-uv-plot', 'value'),
          State('fig-uv-coverage', 'figure'), prevent_initial_call='initial_duplicate')
def update_uv_figure(highlight_antennas: list[str], figure):
    if figure is None:
        raise PreventUpdate

    highlight_colors = [
        "#FF0000",  # Red
        "#0000FF",  # Blue
        "#008000",  # Green
        "#FFA500",  # Orange
        "#800080",  # Purple
        "#00FFFF",  # Cyan
        "#FF00FF",  # Magenta
        "#FFD700",  # Gold
        "#00FF00",  # Lime
        "#A52A2A",  # Brown
        "#FFC0CB",  # Pink
        "#808000",  # Olive
        "#008080",  # Teal
        "#000080",  # Navy
        "#FF7F50",  # Coral
        "#4B0082",  # Indigo
        "#FF8C00",  # DarkOrange
        "#40E0D0",  # Turquoise
        "#6A5ACD",  # SlateBlue
        "#006400",  # DarkGreen
    ]
    if highlight_antennas is not None and len(highlight_antennas) > len(highlight_colors):
        highlight_colors = highlight_colors * (len(highlight_antennas) // len(highlight_colors))

    def get_color(filter_antenna: str):
        if filter_antenna in highlight_antennas:
            return highlight_colors[highlight_antennas.index(filter_antenna)]

        return 'black'

    def get_size(filter_antenna: str):
        if filter_antenna in highlight_antennas:
            return 4

        return 2

    updated_fig = Patch()
    for i, trace in enumerate(figure['data']):
        # trace is a dictionary
        if any(a in highlight_antennas for a in trace['name'].split('-')):
            updated_fig['data'][i]['marker']['color'] = get_color([a for a in trace['name'].split('-')
                                                                   if a in highlight_antennas][0])
            updated_fig['data'][i]['marker']['size'] = 4
        elif updated_fig['data'][i]['marker']['color'] != 'black':
            updated_fig['data'][i]['marker']['color'] = 'black'
            updated_fig['data'][i]['marker']['size'] = 2

    return updated_fig


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

    if duration <= 0.5:
        return 'Must be a positive number', 'form-text text-danger', 'form-control is-invalid'
    elif duration > 4*24:
        return 'Must be shorter than 4 days',  'form-text text-danger', 'form-control is-invalid'

    return "", no_update, 'form-control'
