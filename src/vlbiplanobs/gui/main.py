# from dependencies import Input, Output, State
import os
import random
import threading
from concurrent.futures import ThreadPoolExecutor
from typing import Optional
from datetime import datetime as dt
from dash import Dash, html, dcc, Output, Input, State, no_update
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from astropy import units as u
from astropy.time import Time
from loguru import logger
from vlbiplanobs import sources
from vlbiplanobs import observation
from vlbiplanobs import cli
from vlbiplanobs.gui import inputs, outputs, plots
from vlbiplanobs.gui.callbacks import *
from vlbiplanobs.gui import layout


if os.access("/var/log/planobs.log", os.W_OK):
    logfilename = "/var/log/planobs.log"
else:
    logfilename = "~/log-planobs.log"

_LOG = logger.add(logfilename, backtrace=True, diagnose=True,
                  format="{time:YYYY-MM-DD HH:mm} |  {level} {message}")


class Obs():
    def __init__(self, o: Optional[cli.VLBIObs] = None):
        self.plan_lock: threading.Lock = threading.Lock()
        self._OBS: cli.VLBIObs | None = o
        self.prev_datarate: int = 2048
        self.prev_channels: int = 64
        self.prev_subbands: int = 8

    def set(self, obs: cli.VLBIObs):
        with self.plan_lock:
            self._OBS = obs

    def get(self) -> Optional[cli.VLBIObs]:
        # with self.plan_lock:
        return self._OBS


_main_obs = Obs()

current_directory = os.path.dirname(os.path.realpath(__file__))
external_stylesheets: list = ['https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css']
external_scripts: list = []

app = Dash(__name__, title='EVN Observation Planner', external_scripts=external_scripts,
           external_stylesheets=[dbc.themes.FLATLY, dbc.icons.BOOTSTRAP,
                                 dbc.icons.FONT_AWESOME, dmc.styles.DATES] + external_stylesheets,
           assets_folder=current_directory+'/assets/', eager_loading=False,
           prevent_initial_callbacks=True)


@app.callback([Output('download-data', 'data'),
               Output('downloading', 'children')],
              Input("button-download", "n_clicks"),
              running=[(Output('button-download', 'disabled'), True, False),])
def download_pdf_summary(n_clicks):
    if n_clicks is not None:
        try:
            logger.info("PDF has been requested and created.")
            return dcc.send_bytes(outputs.summary_pdf(_main_obs.get()).getvalue(),
                                  f"planobs_summary-{random.getrandbits(10)}.pdf"), html.Div()
        except ValueError as e:
            print(f"An error occurred: {e}")
            logger.exception("While downloading the PDF: {e}", colorize=True)
            return no_update, no_update

    return no_update, no_update


@app.callback([Output('user-message', 'children'),
               Output('loading-div', 'children'),
               Output('download-summary-div', 'hidden'),
               Output('card-rms', 'children'),
               Output('sensitivity-baseline-modal', 'children'),
               Output('card-resolution', 'children'),
               Output('out-sun', 'children'),
               Output('out-phaseref', 'children'),
               Output('out-ant', 'children'),
               Output('out-elevations', 'hidden'),
               Output('out-elevations-info', 'children'),
               Output('fig-elevations', 'figure'),
               Output('fig-elevations2', 'figure'),
               Output('out-uv-coverage', 'hidden'),
               Output('out-uv-coverage-info', 'children'),
               Output('fig-uv-coverage', 'figure'),
               Output('select-antenna-uv-plot', 'options'),
               Output('div-card-fov', 'children'),
               Output('div-card-vel', 'children'),
               Output('out-worldmap', 'hidden'),
               Output('fig-worldmap', 'figure'),
               Output('card-datasize', 'children')],
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
               State('switch-specify-e-evn', 'value'),
               State('switches-antennas', 'value'),
               [State(f"network-{network}", 'value') for network in observation._NETWORKS]],
              running=[(Output("compute-observation", "disabled"), True, False),],
              suppress_callback_exceptions=True)
def compute_observation(n_clicks, band: int, defined_source: bool, source: str, onsourcetime: float,
                        defined_epoch: bool, startdate: str, duration: float, datarate: int, subbands: int,
                        channels: int, pols: int, inttime: int, e_evn: bool, selected_antennas: list[str],
                        selected_networks: list[bool]):
    """Computes all products to be shown concerning the set observation.
    """
    n_outputs = 22
    if n_clicks is None:
        raise PreventUpdate

    if band == 0 or (not selected_antennas) or (duration is None and (source == '' or not defined_source)):
        return outputs.warning_card("Select the band, antennas, and duration",
                                    "If no source is provided, a duration for the observation "
                                    "must be set."), \
            *[no_update]*(n_outputs - 1)

    selected_antennas = [ant for ant in selected_antennas
                         if observation._STATIONS[ant].has_band(inputs.band_from_index(band))
                         and (not e_evn or observation._STATIONS[ant].real_time)]
    if not selected_antennas:
        return outputs.error_card("No antennas are able to observe with the current setup",
                                  "First, select antennas that can actually observe."), \
               *[no_update]*(n_outputs - 1)

    if defined_epoch and ((startdate is not None and duration is None) or
                          (startdate is None and duration is not None)) == 1:
        return outputs.error_card('The observing epoch is partially defined',
                                  'If you define the observing epoch, then all start date, time, '
                                  'and duration are required.'), *[no_update]*(n_outputs - 1)

    t0 = dt.now()
    network_names = [nn for nb, nn in zip(selected_networks, observation._NETWORKS) if nb]
    try:
        logger.info(f"New Observation: Networks:{','.join(network_names)}; "
                    f"antennas: {','.join(selected_antennas)};"
                    f"band: {inputs.band_from_index(band)}; target: {source}; duration: {duration}h;"
                    f"defined_epoch: {defined_epoch}.")
        _main_obs.set(cli.main(band=inputs.band_from_index(band), stations=sorted(selected_antennas),
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
        logger.exception("An error has occured: {e}.")
        return outputs.error_card("An error has occured", str(e)), *[no_update]*(n_outputs - 1)

    assert _main_obs.get() is not None, "Observation should have been created."
    try:
        futures = {}
        with ThreadPoolExecutor() as executor:
            futures['rms'] = executor.submit(_main_obs.get().thermal_noise)
            futures['beam'] = executor.submit(_main_obs.get().synthesized_beam)
            futures['out-rms'] = executor.submit(outputs.rms, _main_obs.get())
            futures['out-res'] = executor.submit(outputs.resolution, _main_obs.get())
            futures['out_ant'] = executor.submit(outputs.ant_warning, _main_obs.get())
            futures['out_phaseref'] = executor.submit(outputs.warning_low_high_freq,
                                                      _main_obs.get())
            futures['out_fov'] = executor.submit(outputs.field_of_view, _main_obs.get())
            futures['out_freq'] = executor.submit(outputs.summary_freq_res, _main_obs.get())

            out_rms = futures['rms'].result()
            # beam = futures['beam'].result()
            out_rms = futures['out-rms'].result()
            out_res = futures['out-res'].result()
            out_ant = futures['out_ant'].result()
            out_phaseref = futures['out_phaseref'].result()
            out_fov = futures['out_fov'].result()
            out_freq = futures['out_freq'].result()

        with ThreadPoolExecutor() as executor:
            if not (not _main_obs.get().sourcenames or not defined_source):
                # futures['out_plot_elev'] = executor.submit(outputs.plot_elevations, _main_obs.get())
                futures['out_sun'] = executor.submit(outputs.sun_warning, _main_obs.get())
                # futures['out_plot_uv'] = executor.submit(outputsplot_uv_coverage, _main_obs.get())

            # futures['out_worldmap'] = executor.submit(outputsworldmap_plot, _main_obs.get())

            out_sun = no_update if 'out_sun' not in futures else futures['out_sun'].result()
        #     out_plot_uv = ['out_plot_uv' not in futures, no_update if 'out_plot_uv' not in futures \
        #                    else futures['out_plot_uv'].result()]
        #     out_worldmap = futures['out_worldmap'].result()
        #
        #
        # # out_rms = _main_obs.get().thermal_noise()
        # # beam = _main_obs.get().synthesized_beam()
        # # out_rms = outputsrms(_main_obs.get())
        # # out_sens = outputsresolution(_main_obs.get())

        if not _main_obs.get().sourcenames or not defined_source:
            #     out_plot_elev = [True, no_update]
            out_plot_elev = [True, no_update, no_update, no_update]
            #     out_sun = no_update
            out_plot_uv = [True, no_update, no_update, no_update]
        else:
            out_plot_elev = [False,
                             outputs.print_observability_ranges(_main_obs.get()),
                             plots.elevation_plot(_main_obs.get()),
                             plots.elevation_plot_curves(_main_obs.get())]
            #     out_plot_elev = [False, outputsplot_elevations(_main_obs.get())]
            #     out_sun = outputssun_warning(_main_obs.get())
            out_plot_uv = [False, outputs.print_baseline_lengths(_main_obs.get()),
                           plots.uvplot(_main_obs.get()), outputs.put_antenna_options(_main_obs.get())]
        #
        # out_ant = outputsant_warning(_main_obs.get())
        # out_phaseref = outputswarning_phase_referencing_high_freq(_main_obs.get())
        # out_fov = outputsfield_of_view(_main_obs.get())
        # out_freq = outputssummary_freq_res(_main_obs.get())
        out_worldmap = plots.plot_worldmap_stations(_main_obs.get())
        out_baseline_sens = outputs.baseline_sensitivities(_main_obs.get())
        out_datasize = outputs.data_size(_main_obs.get())
    except sources.SourceNotVisible:
        return outputs.error_card('Source Not Visible!',
                                  'The source cannot be observed by the given antennas and/or '
                                  'during the given observing time.'), *[no_update]*(n_outputs - 1)
    except Exception as e:
        logger.exception("While computing: {e}.")
        return outputs.error_card("An error has occured", str(e)), *[no_update]*(n_outputs - 1)

    print(f"Execution time: {(dt.now() - t0).total_seconds()} s")
    logger.info(f"Execution time: {(dt.now() - t0).total_seconds()} s")
    return html.Div(), html.Div(), False, out_rms, out_baseline_sens, out_res, out_sun, \
        out_phaseref, out_ant, *out_plot_elev, *out_plot_uv, out_fov, out_freq, False, out_worldmap, out_datasize


server = app.server
app.index_string = app.index_string.replace('<body>', '<body class="g-sidenav-show bg-gray-100">')

app.layout = dmc.MantineProvider(dbc.Container(fluid=True, className='bg-gray-100 row m-0 p-4', children=[
                   layout.top_banner(app),
                   # inputs.modal_welcome(),
                   html.Div(id='main-window', className='container-fluid d-flex row p-0 m-0',
                            children=[html.Div(id='right-column', className='col-12 col-sm-6 m-0 p-0',
                                               children=layout.inputs_column(app)),
                                      html.Div(id='left-column', className='col-12 col-sm-6 m-0 p-0',
                                               children=[layout.compute_buttons(app),
                                                         layout.outputs_column(app)])]),
                   html.Div(html.A(html.I(className="fa-solid fa-circle-info", style={"font-size": "4rem"}),
                                   id="more-info-button",
                                   className="btn-floating-info btn-lg rounded-circle")),
                   dbc.Tooltip("Opens more information", target='more-info-button'),
                   dbc.Offcanvas(children=inputs.modal_general_info(), id='more-info-modal',
                                 is_open=False, className='shadow-lg blur', placement='end'),
                   html.Div(id='bottom-banner', children=[html.Br(), html.Br(), html.Br()])]))


def main(debug: bool = False):
    return app.run(debug=debug)


if __name__ == '__main__':
    main(debug=True)
