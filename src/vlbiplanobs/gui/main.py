import os
import argparse
from loguru import logger
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime as dt
from dash import Dash, html, dcc, Output, Input, State, no_update, ALL
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from astropy.utils.iers import conf as iers_conf
from astropy import units as u
from astropy.time import Time
from vlbiplanobs import sources
from vlbiplanobs import observation
from vlbiplanobs import cli
from vlbiplanobs.gui import inputs, outputs, plots
from vlbiplanobs.gui.callbacks import *  # noqa: F403
from vlbiplanobs.gui import layout


iers_conf.auto_download = False
iers_conf.auto_max_age = None

if os.access("/var/log/planobs.log", os.W_OK):
    logfilename = "/var/log/planobs.log"
else:
    logfilename = "~/log-planobs.log"

_LOG = logger.add(logfilename, backtrace=True, diagnose=True,
                  format="{time:YYYY-MM-DD HH:mm} |  {level} {message}")

current_directory = os.path.dirname(os.path.realpath(__file__))
external_stylesheets: list = ['https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css']
external_scripts: list = []

app = Dash(__name__, title='EVN Observation Planner', external_scripts=external_scripts,
           external_stylesheets=[dbc.themes.FLATLY, dbc.icons.BOOTSTRAP,
                                 dbc.icons.FONT_AWESOME, dmc.styles.DATES] + external_stylesheets,
           assets_folder=current_directory+'/assets/', eager_loading=False,
           prevent_initial_callbacks=True)


@app.callback([Output('download-data', 'data'),
               Output('downloading', 'children'),
               Output('user-message', 'children', allow_duplicate=True)],
              Input("button-download", "n_clicks"),
              State("store-obs-params", "data"),
              running=[(Output('button-download', 'disabled'), True, False),],
              prevent_initial_call=True)
def download_pdf_summary(n_clicks, obs_params: dict):
    if n_clicks is None or obs_params is None:
        return no_update, no_update, no_update

    try:
        logger.info("PDF generation started on download request.")
        # Recreate observation from stored parameters
        obs = cli.main(
            band=obs_params['band'],
            stations=obs_params['stations'],
            targets=obs_params['targets'],
            duration=obs_params['duration'] * u.h if obs_params['duration'] is not None else None,
            ontarget=obs_params['ontarget'],
            start_time=Time(dt.strptime(f"{obs_params['startdate']} {obs_params['starttime']}", '%Y-%m-%d %H:%M'),
                           format='datetime', scale='utc') if obs_params['startdate'] else None,
            datarate=obs_params['datarate'] * u.Mbit / u.s,
            subbands=obs_params['subbands'],
            channels=obs_params['channels'],
            polarizations=obs_params['polarizations'],
            inttime=obs_params['inttime'] * u.s
        )
        try:
            tmpfile = outputs.summary_pdf(obs, show_figure=True)
        except Exception as fig_error:
            logger.warning(f"Could not include figure in PDF: {fig_error}")
            tmpfile = outputs.summary_pdf(obs, show_figure=False)

        logger.info("PDF has been created.")
        return dcc.send_file(tmpfile, filename="planobs_summary.pdf"), html.Div(), no_update
    except ValueError as e:
        print(f"An error occurred: {e}")
        logger.exception(f"While downloading the PDF: {e}", colorize=True)
        return no_update, no_update, outputs.error_card("Error during the PDF creation",
                                                        "Re-calculate a full observation and try again")


@app.callback([Output('user-message', 'children'),
               Output('loading-div', 'children'),
               Output('download-summary-div', 'hidden'),
               Output('download-link', 'href'),
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
               Output('card-datasize', 'children'),
               Output('div-card-time', 'children'),
               Output('store-prev-datarate', 'data'),
               Output('store-prev-channels', 'data'),
               Output('store-prev-subbands', 'data'),
               Output('store-obs-params', 'data'),
               Output('store-uv-data', 'data')],
              Input('compute-observation', 'n_clicks'),
              [State('band-slider', 'value'),
               State('switch-specify-source', 'value'),
               State('source-input', 'value'),
               State('onsourcetime', 'value'),
               State('switch-specify-epoch', 'value'),
               State('startdate', 'value'),
               State('starttime', 'value'),
               State('duration', 'value'),
               State('datarate', 'value'),
               State('subbands', 'value'),
               State('channels', 'value'),
               State('pols', 'value'),
               State('inttime', 'value'),
               State('switch-specify-e-evn', 'value'),
               State('switches-antennas', 'value'),
               State({'type': 'network-switch', 'index': ALL}, 'value'),
               State({'type': 'network-switch', 'index': ALL}, 'disabled')],
              running=[(Output("compute-observation", "disabled"), True, False),],
              suppress_callback_exceptions=True, prevent_initial_call=True)
@observation.enforce_types
def compute_observation(n_clicks, band: int, defined_source: bool, source: str, onsourcetime: int,
                        defined_epoch: bool, startdate: str, starttime: str, duration: int | float, datarate: int, subbands: int,
                        channels: int, pols: int, inttime: int, e_evn: bool, selected_antennas: list[str],
                        selected_networks: list[bool], disabled_networks: list[bool]):
    """Computes all products to be shown concerning the set observation.
    """
    if n_clicks is None:
        raise PreventUpdate

    vals4error = no_update, True, no_update, *[html.Div()]*6, True, html.Div(), no_update, \
        no_update, True, html.Div(), no_update, no_update, html.Div(), html.Div(), True, no_update, \
        *[html.Div()]*2, no_update, no_update, no_update, no_update, no_update
    if band == 0 or (not selected_antennas) or (duration is None and (source == '' or not defined_source)):
        return outputs.warning_card("Select the band, antennas, and source or duration",
                                    "If no source is provided, a duration for the observation "
                                    "must be set."), *vals4error

    selected_antennas = [ant for ant in selected_antennas
                         if observation._STATIONS[ant].has_band(inputs.band_from_index(band))
                         and (not e_evn or observation._STATIONS[ant].real_time)]
    if not selected_antennas:
        return outputs.error_card("No antennas are able to observe with the current setup",
                                  "First, select antennas that can actually observe."), *vals4error

    if defined_epoch and ((startdate is not None and duration is None) or
                          (startdate is None and duration is not None)):
        return outputs.error_card('The observing epoch is partially defined',
                                  'If you define the observing epoch, then all start date, time, '
                                  'and duration are required.'), *vals4error

    t0 = dt.now()
    network_names = [nn for nb, ns, nn in zip(selected_networks, disabled_networks, observation._NETWORKS) \
                     if nb and not ns]
    try:
        logger.info(f"New Observation: Networks:{','.join(network_names)}; "
                    f"antennas: {','.join(selected_antennas)};"
                    f"datarate: {datarate};",
                    f"band: {inputs.band_from_index(band)}; target: {source}; duration: {duration}"
                    f"{' h' if duration is not None else ''};"
                    f"defined_epoch: {defined_epoch}.")
        obs = cli.main(band=inputs.band_from_index(band), stations=sorted(selected_antennas),
                      targets=[source,] if defined_source and source.strip() != '' else None,
                      duration=duration*u.h if duration is not None else None,
                      ontarget=onsourcetime/100,
                      start_time=Time(dt.strptime(f"{startdate} {starttime}", '%Y-%m-%d %H:%M'),
                                      format='datetime', scale='utc') if defined_epoch else None,
                      datarate=(datarate if not isinstance(datarate, str) else int(datarate))*u.Mbit/u.s,
                      subbands=subbands, channels=channels, polarizations=pols,
                      inttime=inttime*u.s)
        # Pre-compute all cached data BEFORE parallel execution to ensure thread safety.
        # The Observation class methods are now thread-safe with locks, so we can
        # parallelize calls that have no dependencies on each other.
        if obs is not None:
            # Phase 1: Independent methods (no dependencies)
            with ThreadPoolExecutor() as executor:
                f_observable = executor.submit(obs.is_observable)
                f_sun_constraint = executor.submit(obs.sun_constraint)
                f_uv_data = executor.submit(obs.get_uv_data)
                f_baseline_sens = executor.submit(obs.baseline_sensitivity)
                # Wait for phase 1 to complete
                _ = f_observable.result()
                _ = f_sun_constraint.result()
                _ = f_uv_data.result()
                _ = f_baseline_sens.result()

            # Phase 2: Methods that depend on phase 1
            with ThreadPoolExecutor() as executor:
                f_always_obs = executor.submit(obs.is_always_observable)
                f_sun_epochs = executor.submit(obs.sun_limiting_epochs)
                f_thermal = executor.submit(obs.thermal_noise)
                f_uv_values = executor.submit(obs.get_uv_values)
                # Wait for phase 2 to complete
                _ = f_always_obs.result()
                _ = f_sun_epochs.result()
                _ = f_thermal.result()
                _ = f_uv_values.result()

            # Phase 3: Methods that depend on phase 2
            _ = obs.synthesized_beam()
    except ValueError as e:
        logger.exception(f"An error has occured: {e}.")
        return outputs.error_card("Could not plan the observation",
                                  "Your source is not visible during the defined time by >1 antenna."), \
            *vals4error
    except (sources.SourceNotVisible, IndexError):
        return outputs.error_card('Source Not Visible!',
                                  'The source cannot be observed by at least more than one antenna '
                                  'during the given observing time.'), *vals4error
    except Exception as e:
        logger.exception(f"An error has occured: {e}.")
        return outputs.error_card("Could not plan the observation",
                                  "Missing necessary fields in the observation configuration"), *vals4error

    assert obs is not None, "Observation should have been created."
    try:
        if not all(obs.is_observable_by_network(min_stations=1).values()):
            raise sources.SourceNotVisible

        has_source = obs.sourcenames and defined_source
        futures = {}

        # Single ThreadPoolExecutor for all parallel work
        with ThreadPoolExecutor() as executor:
            # Core output functions
            futures['out_rms'] = executor.submit(outputs.rms, obs)
            futures['out_res'] = executor.submit(outputs.resolution, obs)
            futures['out_ant'] = executor.submit(outputs.ant_warning, obs)
            futures['out_phaseref'] = executor.submit(outputs.warning_low_high_freq, obs)
            futures['out_fov'] = executor.submit(outputs.field_of_view, obs)
            futures['out_freq'] = executor.submit(outputs.summary_freq_res, obs)
            futures['out_baseline_sens'] = executor.submit(outputs.baseline_sensitivities, obs)
            futures['out_datasize'] = executor.submit(outputs.data_size, obs)
            futures['out_obstime'] = executor.submit(outputs.obs_time, obs)
            futures['plot_worldmap'] = executor.submit(plots.plot_worldmap_stations, obs)

            # Source-dependent functions
            if has_source:
                futures['out_sun'] = executor.submit(outputs.sun_warning, obs)
                futures['plot_elev'] = executor.submit(plots.elevation_plot, obs)
                futures['plot_elev2'] = executor.submit(plots.elevation_plot_curves, obs)
                futures['plot_uv'] = executor.submit(plots.uvplot, obs)
                futures['uv_data'] = executor.submit(plots.serialize_uv_data, obs)
                futures['out_elev_info'] = executor.submit(outputs.print_observability_ranges, obs)
                futures['out_uv_info'] = executor.submit(outputs.print_baseline_lengths, obs)
                futures['out_ant_options'] = executor.submit(outputs.put_antenna_options, obs)

            # Collect results
            out_rms = futures['out_rms'].result()
            out_res = futures['out_res'].result()
            out_ant = futures['out_ant'].result()
            out_phaseref = futures['out_phaseref'].result()
            out_fov = futures['out_fov'].result()
            out_freq = futures['out_freq'].result()
            out_baseline_sens = futures['out_baseline_sens'].result()
            out_datasize = futures['out_datasize'].result()
            out_obstime = futures['out_obstime'].result()
            out_worldmap = futures['plot_worldmap'].result()

            if has_source:
                out_sun = futures['out_sun'].result()
                out_plot_elev = [False, futures['out_elev_info'].result(),
                                futures['plot_elev'].result(), futures['plot_elev2'].result()]
                out_plot_uv = [False, futures['out_uv_info'].result(),
                            futures['plot_uv'].result(), futures['out_ant_options'].result()]
                uv_data = futures['uv_data'].result()
            else:
                out_sun = no_update
                out_plot_elev = [True, no_update, no_update, no_update]
                out_plot_uv = [True, no_update, no_update, no_update]
                uv_data = no_update
    except sources.SourceNotVisible:
        return outputs.error_card('Source Not Visible!',
                                  'The source cannot be observed by at least more than one antenna '
                                  'during the given observing time.'), *vals4error
    except Exception as e:
        logger.exception(f"While computing: {e}.")
        return outputs.error_card("An error has occured", str(e)), *vals4error

    logger.info(f"Execution time: {(dt.now() - t0).total_seconds()} s")

    # Store observation parameters for lazy PDF generation
    obs_params = {
        'band': inputs.band_from_index(band),
        'stations': sorted(selected_antennas),
        'targets': [source] if defined_source and source.strip() != '' else None,
        'duration': duration,
        'ontarget': onsourcetime / 100,
        'startdate': startdate if defined_epoch else None,
        'starttime': starttime if defined_epoch else None,
        'datarate': datarate if not isinstance(datarate, str) else int(datarate),
        'subbands': subbands,
        'channels': channels,
        'polarizations': pols,
        'inttime': inttime
    }

    return html.Div(), html.Div(), False, None, out_rms, out_baseline_sens, out_res, out_sun, \
        out_phaseref, out_ant, *out_plot_elev, *out_plot_uv, out_fov, out_freq, False, out_worldmap, \
        out_datasize, out_obstime, datarate, channels, subbands, obs_params, uv_data


server = app.server
app.index_string = app.index_string.replace('<body>', '<body class="g-sidenav-show bg-gray-100">')

app.layout = dmc.MantineProvider(dbc.Container(fluid=True, className='bg-gray-100 row m-0 p-4', children=[
                   # Stores for preserving previous values across callbacks (avoids shared state with gunicorn)
                   dcc.Store(id='store-prev-datarate', data=2048),
                   dcc.Store(id='store-prev-channels', data=64),
                   dcc.Store(id='store-prev-subbands', data=8),
                   dcc.Store(id='store-obs-params', data=None),
                   dcc.Store(id='store-uv-data', data=None),
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
    usage = "%(prog)s [-h]  OPTIONS"
    description = "EVN Observation Planner (GUI)"
    parser = argparse.ArgumentParser(description=description, prog="planobs-server", usage=usage)
    parser.add_argument('-d', '--debug', action='store_true', default=False, help="Enable debug mode")
    args = parser.parse_args()
    return app.run(debug=args.debug)


if __name__ == '__main__':
    main(debug=False)
