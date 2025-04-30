# from dash.dependencies import Input, Output, State
import os
import threading
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
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
from vlbiplanobs.gui.callbacks import *
from vlbiplanobs.gui.layout import *


class Obs():
    def __init__(self, o: Optional[cli.VLBIObs] = None):
        self.plan_lock: threading.Lock = threading.Lock()
        self._OBS: cli.VLBIObs | None = o
        self.prev_datarate: int | None = None
        self.prev_channels: int | None = None
        self.prev_subbands: int | None = None

    def set(self, obs: cli.VLBIObs):
        with self.plan_lock:
            self._OBS = obs

    def get(self) -> Optional[cli.VLBIObs]:
        with self.plan_lock:
            return self._OBS


_main_obs = Obs()

current_directory = os.path.dirname(os.path.realpath(__file__))
external_stylesheets: list = []
external_scripts: list = []

app = Dash(__name__, title='EVN Observation Planner', external_scripts=external_scripts,
           external_stylesheets=[dbc.themes.FLATLY, dbc.icons.BOOTSTRAP,
                                 dbc.icons.FONT_AWESOME, dmc.styles.DATES] + external_stylesheets,
           assets_folder=current_directory+'/assets/', eager_loading=True,
           # prevent_initial_callbacks='initial_duplicate')  # type: ignore
           # , suppress_callback_exceptions=True)
           prevent_initial_callbacks=True)



@app.callback([Output('user-message', 'children'),
               Output('loading-div', 'children'),
               Output('button-download-summary', 'children'),
               Output('card-rms', 'children'),
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
               State('switches-antennas', 'value')],
              running=[(Output("compute-observation", "disabled"), True, False),])
def compute_observation(n_clicks, band, defined_source, source, onsourcetime, defined_epoch,
                        startdate, duration, datarate, subbands, channels,
                        pols, inttime, selected_antennas):
    """Computes all products to be shown concerning the set observation.
    """
    if n_clicks is None:
        return [dash.no_update]*16

    if band == 0 or (not selected_antennas) or (duration is None and (source == '' or not defined_source)):
        return out.info_card("Select an observing band and the antennas",
                             "That's the minimum to compute an observation. "
                             "If source is not provided, a duration must be set too."), *[dash.no_update]*15

    if len([ant for ant in selected_antennas if observation._STATIONS[ant]
            .has_band(list(fs.bands.keys())[band-1])]) == 0:
        return out.error_card("No antennas are able to observe at this band",
                              "First, select antennas that can observe at the selected band"), \
               *[dash.no_update]*15

    if defined_epoch and ((startdate is not None and duration is None) or
                          (startdate is None and duration is not None)) == 1:
        return out.error_card('The observing epoch is partially defined',
                              'If you define the observing epoch, all information: start date and time, '
                              'and duration is required.'), *[dash.no_update]*15

    t0 = dt.now()
    try:
        # TODO: Do not create again the Obs. Modify values instead (unless it was None before)
        _main_obs.set(cli.main(band=list(fs.bands.keys())[band-1], stations=selected_antennas,
                      targets=[source,] if defined_source and source.strip() != '' else None,
                      duration=duration*u.h if duration is not None else None,
                      ontarget=onsourcetime/100,
                      start_time=Time(dt.strptime(startdate, '%Y-%m-%d %H:%M:%S'),
                                      format='datetime', scale='utc') if defined_epoch else None,
                      datarate=datarate, subbands=subbands, channels=channels, polarizations=pols,
                      gui=False, tui=False))
        _main_obs.prev_datarate = datarate
        _main_obs.prev_channels = channels
        _main_obs.prev_subbands = subbands
        # I need to run this first otherwise the other functions will fail
        # (likely partially initialized uv values)
        _main_obs.get().thermal_noise()
        _main_obs.get().synthesized_beam()
        _main_obs.get().get_uv_data()
    except Exception as e:
        return out.error_card(f"An error has occured ({e})"), *[dash.no_update]*15

    try:
        futures = {}
        with ThreadPoolExecutor() as executor:
            futures['rms'] = executor.submit(_main_obs.get().thermal_noise)
            futures['beam'] = executor.submit(_main_obs.get().synthesized_beam)
            futures['out-rms'] = executor.submit(out.rms, _main_obs.get())
            futures['out-res'] = executor.submit(out.resolution, _main_obs.get())
            futures['out_ant'] = executor.submit(out.ant_warning, _main_obs.get())
            futures['out_phaseref'] = executor.submit(out.warning_phase_referencing_high_freq,
                                                      _main_obs.get())
            futures['out_fov'] = executor.submit(out.field_of_view, _main_obs.get())
            futures['out_freq'] = executor.submit(out.summary_freq_res, _main_obs.get())

            out_rms = futures['rms'].result()
            beam = futures['beam'].result()
            out_rms = futures['out-rms'].result()
            out_sens = futures['out-res'].result()
            out_ant = futures['out_ant'].result()
            out_phaseref = futures['out_phaseref'].result()
            out_fov = futures['out_fov'].result()
            out_freq = futures['out_freq'].result()

        with ThreadPoolExecutor() as executor:
            if not (not _main_obs.get().sourcenames or not defined_source):
                futures['out_plot_elev'] = executor.submit(out.plot_elevations, _main_obs.get())
                futures['out_sun'] = executor.submit(out.sun_warning, _main_obs.get())
                # futures['out_plot_uv'] = executor.submit(out.plot_uv_coverage, _main_obs.get())

            # futures['out_worldmap'] = executor.submit(out.worldmap_plot, _main_obs.get())

            out_plot_elev = ['out_plot_elev' not in futures, dash.no_update
                             if 'out_plot_elev' not in futures
                             else futures['out_plot_elev'].result()]
            out_sun = dash.no_update if 'out_sun' not in futures else futures['out_sun'].result()
        #     out_plot_uv = ['out_plot_uv' not in futures, dash.no_update if 'out_plot_uv' not in futures \
        #                    else futures['out_plot_uv'].result()]
        #     out_worldmap = futures['out_worldmap'].result()
        #
        #
        # # out_rms = _main_obs.get().thermal_noise()
        # # beam = _main_obs.get().synthesized_beam()
        # # out_rms = out.rms(_main_obs.get())
        # # out_sens = out.resolution(_main_obs.get())

        if not _main_obs.get().sourcenames or not defined_source:
            #     out_plot_elev = [True, dash.no_update]
            #     out_sun = dash.no_update
            out_plot_uv = [True, dash.no_update]
        else:
            #     out_plot_elev = [False, out.plot_elevations(_main_obs.get())]
            #     out_sun = out.sun_warning(_main_obs.get())
            out_plot_uv = [False, out.plot_uv_coverage(_main_obs.get())]
        #
        # out_ant = out.ant_warning(_main_obs.get())
        # out_phaseref = out.warning_phase_referencing_high_freq(_main_obs.get())
        # out_fov = out.field_of_view(_main_obs.get())
        # out_freq = out.summary_freq_res(_main_obs.get())
        out_worldmap = out.worldmap_plot(_main_obs.get())

    except sources.SourceNotVisible:
        return out.error_card('Source Not Visible!',
                              'The source cannot be observed by the given antennas and/or '
                              'during the given observing time.'), *[dash.no_update]*7, False, \
               False, *[dash.no_update]*4
    except Exception as e:
        return out.error_card(f"An error has occured ({e})"), *[dash.no_update]*15

    # print(f"Execution time: {(dt.now() - t0).total_seconds()} s")
    return html.Div(), html.Div(), out.button_summary(_main_obs.get()), out_rms, out_sens, out_sun, \
        out_phaseref, out_ant, *out_plot_elev, False, *out_plot_uv, out_fov, out_freq, out_worldmap


server = app.server
app.index_string = app.index_string.replace(
    '<body>',
    '<body class="g-sidenav-show bg-gray-100">'
)

app.layout = dbc.Container(fluid=True, className='bg-gray-100 row m-0 p-4', children=[
                           el.top_banner(app),
                           el.modal_welcome(),
                           html.Div(id='main-window',
                                    className='container-fluid d-flex row p-0 m-0', children=[
                                       html.Div(id='stations-column',
                                                className='col-lg-6 col-md-6 col-6 p-0 m-0',
                                                children=[html.Div(className='', children=[
                                                    el.card(el.pick_band(fs.bands)),
                                                   *el.networks(app, list(observation._NETWORKS.keys())),
                                                   *el.antenna_list(app, observation._STATIONS)]),
                                                   html.Div(className='col-12',
                                                            children=[html.Div(className='row d-flex',
                                                                               children=html.Div([
                                                                html.Div(children=el.source_selection(),
                                                                         className='col-6',
                                                                         style={'min-width': '250px',
                                                                         'padding-right': '0px'}),
                                                       html.Div(children=el.epoch_selection(),
                                                                className='col-6',
                                                                style={'min-width': '300px',
                                                                       'padding-left': '0px'}),
                                                   ], className='col-12 d-flex',
                                                   style={'align-items': 'stretch'})),
                                                   html.Div(children=el.correlations(),
                                                            className='col-12 col-lg-12 col-md-12 '
                                                                      'col-sm-12 col-xs-12')
                                                   ]),
                                         ]),
                                html.Div(id='left-column', className='col-lg-6 col-md-6 col-6 m-0 p-0',
                                         children=[html.Div(className='container row flex', children=[
                                            el.compute_button(),
                                            html.Div(className='col-6 text-center',
                                                     id='button-download-summary',
                                                     children=html.Div(id='button-download'))]),
                                         dcc.Download(id="download-data"),
                                         html.Div(className='col-12 mx-0 px-4 py-2', id='user-message',
                                                  children=[out.info_card("Set your VLBI observation "
                                                      "and press 'calculate'", "The expected outcome "
                                                      "details will appear here.")]),
                                         html.Div(className='col-12 mx-0 pb-2', children=[
                                             html.Div(className='row d-flex m-0', children=[
                                                 html.Div(className='col-6 m-0 px-2 d-flex', id='card-rms',
                                                          children=out.rms()),
                                                 html.Div(className='col-6 m-0 px-2 d-flex', id='card-resolution',
                                                          children=out.resolution())
                                              ], style={'align-items': 'stretch'})
                                         ]),
                                         html.Div(className='col-12 mx-0 mb-2', id='out-sun', children=[]),
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
                                 html.Div(html.A(html.I(className="fa-solid fa-circle-info",
                                                 style={"font-size": "4rem"}),
                                          id="more-info-button",
                                          className="btn-floating-info btn-lg rounded-circle")),
                                 dbc.Tooltip("Opens more information", target='more-info-button'),
                                 dbc.Offcanvas(children=el.modal_general_info(), id='more-info-modal',
                                               is_open=False, className='shadow-lg blur', placement='end')
                           # dbc.Modal([
                           #      dbc.ModalHeader("More Information"),
                           #      dbc.ModalBody("This is some additional information"),
                           #      dbc.ModalFooter(
                           #          dbc.Button("Close", id="close-info-modal", className="ml-auto")
                           #      ),
                           #  ], id="more-info-modal", is_open=False),
                           ]),
                           html.Div(id='bottom-banner', children=[html.Br(), html.Br(), html.Br()])])


def main(debug: bool = False):
    return app.run(debug=debug)


if __name__ == '__main__':
    main(debug=True)
