# from dash.dependencies import Input, Output, State
import os
import threading
from dataclasses import dataclass
from typing import Optional, Union
from datetime import datetime as dt
import numpy as np
import dash
from dash import Dash, html, dcc, callback, Output, Input, State
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from astropy import coordinates as coord
from astropy import units as u
from astropy.time import Time
from vlbiplanobs import freqsetups as fs
from vlbiplanobs import stations
from vlbiplanobs import sources
from vlbiplanobs import observation
from vlbiplanobs import cli
from vlbiplanobs.gui import inputs as el
from vlbiplanobs.gui import outputs as out
from vlbiplanobs.gui import plots
# adding the possibility of disabled. Will be implemented in a future version of dash_bootstrap_components


class Obs():
    _OBS: Optional[cli.VLBIObs] = None

    def __init__(self, o: Optional[cli.VLBIObs] = None):
        self.plan_lock = threading.Lock()
        self._OBS = o
        self.prev_datarate = None

    def set(self, obs: cli.VLBIObs):
        with self.plan_lock:
            self._OBS = obs

    def get(self) -> Optional[cli.VLBIObs]:
        with self.plan_lock:
            return self._OBS

_main_obs = Obs()

current_directory = os.path.dirname(os.path.realpath(__file__))
# external_stylesheets: list = [current_directory + '/assets/css/style-planobs-v3.css']
external_stylesheets: list = []
# current_directory + '/assets/css/style-planobs-v3.css']
# external_stylesheets: list = [current_directory + '/assets/css/' + css \
#                               for css in ('nucleo-icons.css', 'nucleo-svg.css', 'soft-ui-dashboard.css',
#                                           'soft-ui-dashboard.min.css', 'soft-ui-dashboard.css.map')]
external_scripts: list = []
# "https://kit.fontawesome.com/69c65a0ab5.js"]  # [current_directory + '/assets/js/' + js \
# external_scripts: list = [current_directory + '/assets/js/' + js \
#                               for js in ('soft-ui-dashboard.js', 'soft-ui-dashboard.min.js',
#                                          'soft-ui-dashboard.js.map')]  # ,
# "https://planobs.jive.eu/assets/mixpanel-analytics.js"]
for css in external_stylesheets:
    assert os.path.isfile(css), f"The file {css} does not exist"

for js in external_scripts:
    assert os.path.isfile(js), f"The file {js} does not exist"

app = Dash(__name__, title='EVN Observation Planner', external_scripts=external_scripts,
           external_stylesheets=[dbc.themes.FLATLY, dbc.icons.BOOTSTRAP,
                                 dbc.icons.FONT_AWESOME, dmc.styles.DATES] + external_stylesheets,
           assets_folder=current_directory+'/assets/', eager_loading=True,
           prevent_initial_callbacks='initial_duplicate')  # , suppress_callback_exceptions=True)


@app.callback([Output('band-slider', 'marks'),
               *[Output(f"network-{network}-label-wav", 'hidden') for network in observation._NETWORKS],
               *[Output(f"network-{network}-label-freq", 'hidden') for network in observation._NETWORKS]],
              Input('switch-band-label', 'value'))
def change_band_labels(show_wavelengths: bool):
    top_labels = [{'label': 'ν (GHz)', 'style': {'font-weight': 'bold', 'color': '#004990'}}] + \
                 [b.split('or')[1].replace('GHz', '').strip() for b in fs.bands.values()]
    bottom_labels = [{'label': 'λ (cm)', 'style': {'font-weight': 'bold', 'color': '#004990'}}] + \
                    [b.split('or')[0].replace('cm', '').strip() for b in fs.bands.values()]
    labels = {i: l for i, l in enumerate(bottom_labels)} if show_wavelengths \
              else {i: l for i, l in enumerate(top_labels)}
    return labels, \
           *[not show_wavelengths for _ in observation._NETWORKS], \
           *[show_wavelengths for _ in observation._NETWORKS]


@app.callback([Output(f"network-{network}", 'disabled') for network in observation._NETWORKS] + \
              [Output(f"network-{network}-card", 'style') for network in observation._NETWORKS],
              Input('band-slider', 'value'),
              [State(f"network-{network}-card", 'style') for network in observation._NETWORKS])
def enable_networks_with_band(band_index: int, *card_styles):
    if band_index == 0:
        card_styles = [{k: v if k != 'opacity' else 1.0 for k, v in card_style.items()}
                       for card_style in card_styles]
        return [False for _ in observation._NETWORKS] + card_styles

    opacity = lambda x: 1.0 if x else 0.2
    card_styles = [{k: v if k != 'opacity' else \
                    opacity(list(fs.bands.keys())[band_index - 1] in \
                    [n for n in observation._NETWORKS.values()][i].observing_bands) \
                    for k, v in card_style.items()} for i, card_style in enumerate(card_styles)]
    return [list(fs.bands.keys())[band_index - 1] not in observation._NETWORKS[network].observing_bands \
            for network in observation._NETWORKS] + card_styles


@app.callback(Output('datarate', 'value'),
              [Input('band-slider', 'value'),
               *[Input(f"network-{network}", 'value') for network in observation._NETWORKS]])
def update_datarate(band, *networks):
    if not [n for n in networks if n]:
        return dash.no_update

    if band == 0:
        return dash.no_update

    if 'EVN' in networks:
        return cli._NETWORKS['EVN'].max_datarate(list(fs.bands.keys())[band-1])

    try:
        return min([nn.max_datarate(list(fs.bands.keys())[band-1]).value \
                    for ni, nn in zip(networks, cli._NETWORKS.values()) \
                    if ni and list(fs.bands.keys())[band-1] in nn.observing_bands])
    except (KeyError, ValueError):
        return dash.no_update
    return f"{val.value:.01f} {val.unit.to_string("unicode")}"


@app.callback(Output(f"switches-antennas", 'children'),
              [Input('band-slider', 'value'),
               Input('switch-specify-e-evn', 'value')],
              State(f"switches-antennas", 'value'))
def enable_antennas_with_band(band_index: int, do_e_evn: bool, selection_antennas: list):
    if band_index == 0:
        return [el.antenna_card_hover(app, dmc.Chip(s.name, value=s.codename,
                                                    color='#004990', styles={'display': 'grid',
                                             'grid-template-columns': 'repeat(auto-fit, minmax(10rem, 1fr))'},), s) for s in observation._STATIONS]
        return [dmc.Chip(s.name, value=s.codename, color='#004990',
                         styles={'display': 'grid',
                                 'grid-template-columns': 'repeat(auto-fit, minmax(10rem, 1fr))'}) \
                for s in observation._STATIONS] #antennas]

    return [el.antenna_card_hover(app,
                                  dmc.Chip(s.name, value=s.codename,
                                           color='#004990', styles={'display': 'grid',
                                             'grid-template-columns': 'repeat(auto-fit, minmax(10rem, 1fr))'},
                                           disabled=not list(fs.bands.keys())[band_index - 1] in s.bands or \
                                                    (do_e_evn and not s.real_time),
                                           ), s) for s in observation._STATIONS]


@app.callback(Output('switches-antennas', 'value'),
              [Input(f"network-{network}", 'value') for network in observation._NETWORKS],
              State('switches-antennas', 'value'))
def update_selected_antennas_from_networks(*args):
    networks, current_antennas = args[:-1], args[-1]
    ants2exclude, ants2include = [], []
    for ni, nn in zip(networks, cli._NETWORKS.values()):
        if ni:
            ants2include += nn.station_codenames
        else:
            ants2exclude += nn.station_codenames

    ants2exclude = [ant for ant in ants2exclude if ant not in ants2include]
    current_antennas = [ant for ant in current_antennas if ant not in ants2exclude]
    return current_antennas + [ant for ant in ants2include if ant not in current_antennas]


@app.callback(Output('source-selection-div', 'hidden'), Input('switch-specify-source', 'value'))
def toggle_source_field(pick_source: bool):
    return not pick_source


@app.callback(Output('epoch-selection-div', 'hidden'), Input('switch-specify-epoch', 'value'))
def toggle_epoch_field(define_epoch: bool):
    return not define_epoch


@app.callback(Output('accordion-ant', 'className'), Input('accordion-state', 'active_item'))
def toggle_accordion_arrow(accordion_expanded: int):
    return 'fa fa-solid fa-angle-down' if accordion_expanded is None else 'fa fa-solid fa-angle-up'


@app.callback([Output('datarate', 'value', allow_duplicate=True),
               Output('channels', 'value'),
               Output('subbands', 'value')],
              Input('switch-specify-continuum', 'value'),
              [State('band-slider', 'value'),
               State('datarate', 'value'),
               *[State(f"network-{network}", 'value') for network in observation._NETWORKS]])
def prioritize_spectral_line(do_spectral_line: bool, band: str, datarate: int, *networks):
    if band == 0:
        dt = 2048

    if 'EVN' in networks:
        dt = cli._NETWORKS['EVN'].max_datarate(list(fs.bands.keys())[band-1])
    elif not networks:
        dt = min([nn.max_datarate(list(fs.bands.keys())[band-1]).value \
                  for ni, nn in zip(networks, cli._NETWORKS.values()) if ni])
    else:
        dt = 2048

    if do_spectral_line:
        _main_obs.prev_datarate = datarate
    else:
        _main_obs.prev_datarate = dt

    return 32 if do_spectral_line else _main_obs.prev_datarate, 4096 if do_spectral_line else 64, \
           1 if do_spectral_line else 8


@app.callback([Output('error_source', 'children'),
        Output('error_source', 'className'), Output('source-input', 'className')],
        [Input('source-input', 'value')])
def get_initial_source(source_coord):
    """Verifies that the introduced source coordinates have a right format.
    If they are correct, it does nothing. If they are incorrect, it shows an error label.
    """

    if source_coord != 'hh:mm:ss dd:mm:ss' and source_coord != None and source_coord != '':
        if len(source_coord) > 40:
            # Otherwise the source name check gets too slow
            return "Name too long.", 'form-text text-danger', 'form-control'
        try:
            src = sources.Source.source_from_str(source_coord)
            return src.coord.to_string('hmsdms', precision=3), 'form-text', 'form-control is-valid'
        except ValueError as e:
            return "Wrong coordinates.", 'form-text text-danger', 'form-control is-invalid'
        except coord.name_resolve.NameResolveError as e:
            return "Unrecognized name. Use 'hh:mm:ss dd:mm:ss' or 'XXhXXmXXs XXdXXmXXs'", \
                   'form-text text-danger', 'form-control is-invalid'
    else:
        return '', dash.no_update, 'form-control'


@app.callback(Output('bandwidth-label', 'children'),
              [Input('datarate', 'value'), Input('pols', 'value')])
def update_bandwidth_label(datarate: int, npols: int):
    """Updates the total bandwidth label as a function of the selected datarate and number of
    polarizations. Returns a string with the value and units.
    """
    if (None not in (datarate, npols)) and (datarate != -1):
        # Either 1 or 2 pols per station:
        temp = npols % 3 + npols // 3
        return [f"The maximum bandwidth is {cli.optimal_units(datarate*u.MHz/(temp*2*2),
                                                             [u.GHz, u.MHz, u.kHz] )}."
               ]

    return ''


@app.callback([Output("sensitivity-baseline-modal", "is_open"),
               Output("sensitivity-baseline-modal", "children"),
               Output("sens-baseline-style", "style")],
              Input("button-sensitivity-baseline", "n_clicks"),
              State("sensitivity-baseline-modal", "is_open"))
def toggle_modal(n_clicks, is_open):
    if n_clicks is None:
        return False, dash.no_update, {'display': 'none'}

    return not is_open, out.show_baseline_sensitivities(_main_obs.get()), {'display': 'block'}


@app.callback(Output('fig-coverage', 'figure'),
              Input('select-antenna-uv-plot', 'value'))
def update_uv_figure(highlight_antennas: list[str]):
    return plots.uvplot(_main_obs.get(), highlight_antennas)


@app.callback(Output('error_duration', 'children'),
              Input('duration', 'value'))
def check_initial_obstime(duration):
    """Verify the introduced times/dates for correct values.
    Once the user has introduced all values for the start and end of the observation,
    it guarantees that they have the correct shape:
        - the duration of the observation is > 0 hours.
        - The total observing length is less than five days (value chosen for computational reasons).
    """
    if duration is None:
        return ""

    if (not isinstance(duration, float)) and (not isinstance(duration, int)):
        return "Must be a number"

    if duration <= 0.0:
        return "The duration must be a positive number"
    elif duration > 4*24:
        return "Please, put an observation shorter than 4 days"

    return ""


@app.callback([Output('user-message', 'children'),
               Output('card-rms', 'hidden'),
               Output('rms-value', 'children'),
               Output('rms-per-channel-value', 'children'),
               Output('rms-per-time-value', 'children'),
               Output('table-sensitivities', 'children'),
               # Output('table-sensitivities', 'hidden'),  # TODO: maybe this one for hidden?
               Output('card-resolution', 'hidden'),
               Output('card-resolution', 'children'),
               Output('out-sun', 'children'),
               Output('out-phaseref', 'children'),
               Output('out-ant', 'children'),
               Output('out-elevations', 'hidden'),
               Output('out-elevations', 'children'),
               Output("sensitivity-baseline-modal", "is_open", allow_duplicate=True),
               Output('out-uv-coverage', 'hidden'),
               Output('out-uv-coverage', 'children'),
               Output('div-card-fov', 'children'),
               Output('div-card-vel', 'children'),
               Output('out-worldmap', 'children')],
              Input('compute-observation', 'n_clicks'),
              [State('band-slider', 'value'),
               State('switch-specify-source', 'value'),
               State('source-input', 'value'),
               State('onsourcetime', 'value'),
               State('switch-specify-epoch', 'value'),
               State('starttime', 'value'),
               State('duration', 'value'),
               State('datarate', 'value'),
               State('subbands', 'value'),
               State('channels', 'value'),
               State('pols', 'value'),
               State('inttime', 'value'),
               State('switches-antennas', 'value')])
def compute_observation(n_clicks, band, defined_source, source, onsourcetime, defined_epoch,
                        startdate, duration, datarate, subbands, channels,
                        pols, inttime, selected_antennas):
    """Computes all products to be shown concerning the set observation.
    """
    if n_clicks is None:
        return *[dash.no_update]*13, False, *[dash.no_update]*5

    if band == 0 or not selected_antennas or (duration is None and (source == '' or not defined_source)):
        return out.info_card("Select an observing band and the antennas",
                             "That's the minimum to compute an observation"), True, *[dash.no_update]*4, True, \
               *[dash.no_update]*4, True, dash.no_update, False, *[dash.no_update]*5

    if defined_epoch and ((startdate is not None and duration is None) or \
                          (startdate is None and duration is not None)) == 1:
        return out.error_card('The observing epoch is partially defined',
                              'If you define the observing epoch, all information: start date and time, '
                              'and duration is required.'), True, *[dash.no_update]*4, True, \
               *[dash.no_update]*4, True, dash.no_update, False, *[dash.no_update]*5


    obs = cli.main(band=list(fs.bands.keys())[band-1], stations=selected_antennas,
                   targets=[source,] if defined_source and source.strip() != '' else None,
                   duration=duration*u.h if duration is not None else None,
                   ontarget=onsourcetime/100,
                   start_time=Time(dt.strptime(startdate, '%Y-%m-%d %H:%M:%S'),
                                   format='datetime', scale='utc') if defined_epoch else None,
                   datarate=datarate, subbands=subbands, channels=channels, polarizations=pols,
                   gui=False, tui=False)
    _main_obs.set(obs)

    if _main_obs.get().datarate is not None: # and _main_obs.get().duration is not None:
        if isinstance(thermal_noise := _main_obs.get().thermal_noise(), dict):
            rms = list(thermal_noise.values())[0]
        else:
            rms = thermal_noise

        out_rms = [False, out.quantity2str(cli.optimal_units(rms,
                                                             [u.MJy/u.beam, u.kJy/u.beam, u.Jy/u.beam,
                                                              u.mJy/u.beam, u.uJy/u.beam])),
                   out.quantity2str(cli.optimal_units(rms/np.sqrt(1*u.min/ \
                                                      (_main_obs.get().duration \
                                                      if _main_obs.get().duration is not None else 24*u.h)),
                                                      [u.MJy/u.beam, u.kJy/u.beam, u.Jy/u.beam,
                                                       u.mJy/u.beam, u.uJy/u.beam])),
                   out.quantity2str(cli.optimal_units(rms*np.sqrt(_main_obs.get().subbands* \
                                                      _main_obs.get().channels),
                                                      [u.MJy/u.beam, u.kJy/u.beam, u.Jy/u.beam,
                                                       u.mJy/u.beam, u.uJy/u.beam])),
                   out.show_baseline_sensitivities(_main_obs.get())]
    else:
        out_rms = [True, dash.no_update, dash.no_update, dash.no_update, dash.no_update]

    try:
        beam = list(_main_obs.get().synthesized_beam().values())[0]
        sens_units = cli.optimal_units(beam['bmaj'], [u.mas, u.arcsec, u.arcmin, u.deg]).unit
        out_sens = [False, f"{cli.optimal_units(beam['bmaj'], [u.mas, u.arcsec, u.arcmin, u.deg]).value} x " \
                    f"{beam['bmin'].to(sens_units)}2, PA = {beam['pa']}."]
        out_sens = [False, out.resolution(_main_obs.get())]

        if not _main_obs.get().sourcenames or not defined_source:
            out_plot_elev = [True, dash.no_update]
            out_sun = dash.no_update
            out_plot_uv = [True, dash.no_update]
        else:
            out_plot_elev = [False, out.plot_elevations(_main_obs.get())]
            out_sun = out.sun_warning(_main_obs.get())
            out_plot_uv = [False, out.plot_uv_coverage(_main_obs.get())]

        out_ant = out.ant_warning(_main_obs.get())
        out_phaseref = out.warning_phase_referencing_high_freq(_main_obs.get())
        out_fov = out.field_of_view(_main_obs.get())
        out_freq = out.summary_freq_res(_main_obs.get())
    except sources.SourceNotVisible:
        return out.error_card('Source Not Visible!',
                              'The source cannot be observed by the given antennas and/or '
                              'during the given observing time.'), True, *[dash.no_update]*4, True, \
               *[dash.no_update]*4, True, dash.no_update, False, *[dash.no_update]*5

    return html.Div(), *out_rms, *out_sens, out_sun, out_phaseref, out_ant, *out_plot_elev, False, \
           *out_plot_uv, out_fov, out_freq, out.worldmap_plot(_main_obs.get())


# app.config.suppress_callback_exceptions = True  # Avoids error messages for id's that haven't been loaded yet
server = app.server
app.index_string = app.index_string.replace(
    '<body>',
    '<body class="g-sidenav-show bg-gray-100">'
)

app.layout = dbc.Container(fluid=True, className='bg-gray-100 row m-0 p-4',
                           children=el.top_banner(app) + [
                           html.Div(id='main-window', className='container-fluid d-flex row p-0 m-0',
                                    children=[
                                html.Div(id='stations-column', className='col-lg-6 col-md-6 col-6 p-0 m-0',
                                         children=[html.Div(className='',
                                                            children=[el.card(el.pick_band(fs.bands)),
                                                   *el.networks(app, list(observation._NETWORKS.keys()),
                                                                bands=fs.bands),
                                                   *el.antenna_list(app, observation._STATIONS)]),
                                                   html.Div(className='col-12',
                                                            children=[html.Div(className='row',
                                                            children=[
                                                       html.Div(children=el.source_selection(),
                                                                className='col-6 col-lg-6 col-md-6 '
                                                                          'col-sm-6 col-xs-12',
                                                                style={'min-width': '250px',
                                                                       'padding-right': '0px'}),
                                                       html.Div(children=el.epoch_selection(),
                                                                className='col-6 col-lg-6 col-md-6 '
                                                                          'col-sm-6 col-xs-12',
                                                                style={'min-width': '300px',
                                                                       'padding-left': '0px'}),
                                                   ]),
                                                   html.Div(children=el.correlations(),
                                                            className='col-12 col-lg-12 col-md-12 '
                                                                      'col-sm-12 col-xs-12')
                                                   ]),
                                         ]),
                                html.Div(id='left-column', className='col-lg-6 col-md-6 col-6 m-0 p-0',
                                         children=[html.Div(className='row',
                                                            children=el.compute_button()),
                    # dcc.Loading(id="loading2", children=[html.Div(id="loading-output2"), html.Br()],
                                # type="dot"),
                                         html.Div(className='col-12 mx-0', id='user-message',
                                                  children=[out.info_card("Specify your VLBI observation "
                                                      "and then press 'compute'", "The expected outcome "
                                                      "details will appear here.")]),
                                         html.Div(className='row', children=el.results()),
                                         html.Div(className='col-12 mx-0', children=[
                                             html.Div(className='row d-flex m-0', children=[
                                                 html.Div(className='col-6 m-0 px-2', id='card-rms',
                                                          children=out.rms(), hidden=True),
                                                 html.Div(className='col-6 m-0 px-2', id='card-resolution',
                                                          children=out.resolution(), hidden=True)
                                              ], style={'align-items': 'stretch'})
                                         ]),
                                         html.Div(className='col-12 mx-0', id='out-sun', children=[]),
                                         html.Div(className='col-12 mx-0', id='out-phaseref', children=[]),
                                         html.Div(className='col-12 mx-0', id='out-ant', children=[]),
                                         html.Div(className='col-12 mx-0', id='out-elevations',
                                                  children=[]),
                                         html.Div(className='col-12 mx-0', id='out-uv-coverage',
                                                  children=[dcc.Graph(id='fig-coverage'),
                                                            dbc.Select(id="select-antenna-uv-plot")],
                                                  hidden=True),
                                         html.Div(className='row', children=[
                                             html.Div(className='col-12 mx-0', children=[
                                                 html.Div(className='row d-flex m-0', children=[
                                             html.Div(className='col-6 m-0 px-2', id='div-card-fov',
                                                      children=[]),
                                             html.Div(className='col-6 m-0 px-2', id='div-card-vel',
                                                      children=[]),
                                         ])])]),
                                         html.Div(className='row', children=[
                                             html.Div(className='col-12 mx-0', children=[
                                                 html.Div(className='row col-12 mx-0 my-2 px-2',
                                                          id='out-worldmap', children=[]),
                                         ])])
                                 ]),
                           ]),
                           html.Div(id='bottom-banner', children=[html.Br(), html.Br(), html.Br()])])


def main(debug: bool = False):
    return app.run(debug=debug)


if __name__ == '__main__':
    main(debug=True)
