import os
import argparse
from typing import Optional
import json
import importlib.metadata
from packaging.version import Version
from urllib.parse import quote, unquote
from loguru import logger
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime as dt
from dash import Dash, html, dcc, Output, Input, State, no_update, ALL, \
    clientside_callback
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from astropy.utils.iers import conf as iers_conf
from astropy import units as u
from astropy.time import Time
from furl import furl
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
    logfilename = os.path.expanduser("~/log-planobs.log")

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
            start_time=Time(dt.strptime(f"{obs_params['startdate'][:10]} {obs_params['starttime']}", '%Y-%m-%d %H:%M'),
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
               Output('out-uv-coverage', 'hidden'),
               Output('out-uv-coverage-info', 'children'),
               Output('select-antenna-uv-plot', 'options'),
               Output('div-card-fov', 'children'),
               Output('div-card-vel', 'children'),
               Output('out-worldmap', 'hidden'),
               Output('card-datasize', 'children'),
               Output('div-card-time', 'children'),
               Output('store-prev-datarate', 'data'),
               Output('store-prev-channels', 'data'),
               Output('store-prev-subbands', 'data'),
               Output('store-obs-params', 'data'),
               Output('store-uv-data', 'data'),
               Output('store-elev-data', 'data'),
               Output('store-worldmap-data', 'data')],
              Input('compute-observation', 'n_clicks'),
              [State('band-slider', 'value'),
               State('switch-specify-source', 'value'),
               State('source-input', 'value'),
               State('onsourcetime', 'value'),
               State('switch-specify-epoch', 'value'),
               State('startdate', 'date'),
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
def compute_observation(n_clicks, band: int, defined_source: bool, source: Optional[str],
                        onsourcetime: Optional[int], defined_epoch: bool, startdate: Optional[str],
                        starttime: Optional[str], duration: Optional[int | float],
                        datarate: Optional[int], subbands: Optional[int], channels: Optional[int],
                        pols: Optional[int], inttime: Optional[int], e_evn: bool,
                        selected_antennas: list[str], selected_networks: list[bool],
                        disabled_networks: list[bool]):
    """Computes observation and returns text outputs immediately.

    Plot figures are rendered progressively via store-triggered callbacks
    (store-elev-data, store-uv-data, store-worldmap-data) so the user sees
    text results without waiting for heavy Plotly figure generation.
    """
    if n_clicks is None:
        raise PreventUpdate

    # Error return: 27 outputs total (text + stores, no figures)
    vals4error = no_update, True, no_update, *[html.Div()]*6, True, html.Div(), True, html.Div(), \
        no_update, html.Div(), html.Div(), True, *[html.Div()]*2, \
        no_update, no_update, no_update, no_update, no_update, no_update, no_update
    if band == 0 or (not selected_antennas) or (duration is None and (not source or not defined_source)):
        return outputs.warning_card("Select the band, antennas, and source or duration",
                                    "If no source is provided, a duration for the observation "
                                    "must be set."), *vals4error

    if any(v is None for v in (datarate, subbands, channels, pols, inttime, onsourcetime)):
        return outputs.error_card("Missing observation parameters",
                                  "Please set the data rate, subbands, channels, polarizations, "
                                  "integration time, and on-source percentage."), *vals4error

    assert datarate is not None and subbands is not None and channels is not None
    assert pols is not None and inttime is not None and onsourcetime is not None

    selected_antennas = [ant for ant in selected_antennas
                         if observation._STATIONS[ant].has_band(inputs.band_from_index(band))
                         and (not e_evn or observation._STATIONS[ant].real_time)]
    if not selected_antennas:
        return outputs.error_card("No antennas are able to observe with the current setup",
                                  "First, select antennas that can actually observe."), *vals4error

    if duration is not None and duration < 0.1:
        return outputs.error_card('Duration too short',
                                  'Minimum duration is 0.1 h, but you can check the instantaneous '
                                  'sensitivity in the provided values (they are also calculated per '
                                  'minute integration).'), *vals4error

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
                      targets=[source] if defined_source and source and source.strip() != '' else None,
                      duration=duration*u.h if duration is not None else None,
                      ontarget=onsourcetime/100,
                      start_time=Time(dt.strptime(f"{startdate[:10]} {starttime}", '%Y-%m-%d %H:%M'),
                                      format='datetime', scale='utc') if defined_epoch else None,
                      datarate=(datarate if not isinstance(datarate, str) else int(datarate))*u.Mbit/u.s,
                      subbands=subbands, channels=channels, polarizations=pols,
                      inttime=inttime*u.s)
        # Pre-compute all cached data BEFORE parallel execution to ensure thread safety.
        if obs is not None:
            with ThreadPoolExecutor() as executor:
                # Submit all independent pre-computations in a single pool
                pre_futures = [executor.submit(obs.is_observable),
                               executor.submit(obs.sun_constraint),
                               executor.submit(obs.get_uv_data),
                               executor.submit(obs.baseline_sensitivity),
                ]
                for f in pre_futures:
                    f.result()
                # Second wave depends on is_observable being cached
                dep_futures = [executor.submit(obs.is_always_observable),
                               executor.submit(obs.sun_limiting_epochs),
                               executor.submit(obs.thermal_noise),
                               executor.submit(obs.get_uv_values),
                               executor.submit(obs.synthesized_beam),
                ]
                for f in dep_futures:
                    f.result()
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

        with ThreadPoolExecutor() as executor:
            # Text output functions (fast)
            futures['out_rms'] = executor.submit(outputs.rms, obs)
            futures['out_res'] = executor.submit(outputs.resolution, obs)
            futures['out_ant'] = executor.submit(outputs.ant_warning, obs)
            futures['out_phaseref'] = executor.submit(outputs.warning_low_high_freq, obs)
            futures['out_fov'] = executor.submit(outputs.field_of_view, obs)
            futures['out_freq'] = executor.submit(outputs.summary_freq_res, obs)
            futures['out_baseline_sens'] = executor.submit(outputs.baseline_sensitivities, obs)
            futures['out_datasize'] = executor.submit(outputs.data_size, obs)
            futures['out_obstime'] = executor.submit(outputs.obs_time, obs)

            # Lightweight data serialization for deferred plot rendering
            futures['worldmap_data'] = executor.submit(plots.serialize_worldmap_data, obs)  # type: ignore[arg-type]

            if has_source:
                futures['out_sun'] = executor.submit(outputs.sun_warning, obs)
                futures['elev_data'] = executor.submit(plots.serialize_elevation_data, obs)  # type: ignore[arg-type]
                futures['uv_data'] = executor.submit(plots.serialize_uv_data, obs)  # type: ignore[arg-type]
                futures['out_elev_info'] = executor.submit(outputs.print_observability_ranges, obs)
                futures['out_uv_info'] = executor.submit(outputs.print_baseline_lengths, obs)
                futures['out_ant_options'] = executor.submit(outputs.put_antenna_options, obs)  # type: ignore[arg-type]

            # Collect text results
            out_rms = futures['out_rms'].result()
            out_res = futures['out_res'].result()
            out_ant = futures['out_ant'].result()
            out_phaseref = futures['out_phaseref'].result()
            out_fov = futures['out_fov'].result()
            out_freq = futures['out_freq'].result()
            out_baseline_sens = futures['out_baseline_sens'].result()
            out_datasize = futures['out_datasize'].result()
            out_obstime = futures['out_obstime'].result()
            worldmap_data = futures['worldmap_data'].result()

            if has_source:
                out_sun = futures['out_sun'].result()
                out_elev = [False, futures['out_elev_info'].result()]
                out_uv = [False, futures['out_uv_info'].result(), futures['out_ant_options'].result()]
                uv_data = futures['uv_data'].result()
                elev_data = futures['elev_data'].result()
            else:
                out_sun = no_update  # type: ignore[assignment]
                out_elev = [True, no_update]
                out_uv = [True, no_update, no_update]
                uv_data = no_update  # type: ignore[assignment]
                elev_data = no_update  # type: ignore[assignment]
    except sources.SourceNotVisible:
        return outputs.error_card('Source Not Visible!',
                                  'The source cannot be observed by at least more than one antenna '
                                  'during the given observing time.'), *vals4error
    except Exception as e:
        logger.exception(f"While computing: {e}.")
        return outputs.error_card("An error has occured", str(e)), *vals4error

    logger.info(f"Execution time: {(dt.now() - t0).total_seconds()} s")

    obs_params = {
        'band': inputs.band_from_index(band),
        'stations': sorted(selected_antennas),
        'targets': [source] if defined_source and source and source.strip() != '' else None,
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
        out_phaseref, out_ant, *out_elev, *out_uv, out_fov, out_freq, False, \
        out_datasize, out_obstime, datarate, channels, subbands, obs_params, \
        uv_data, elev_data, worldmap_data


server = app.server
app.index_string = app.index_string.replace('<head>', '<head><meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">')
app.index_string = app.index_string.replace('<body>', '<body class="g-sidenav-show bg-gray-100">')

app.layout = dmc.MantineProvider(dbc.Container(fluid=True, className='bg-gray-100 row m-0 p-4', children=[
                   # Stores for preserving previous values across callbacks (avoids shared state with gunicorn)
                   dcc.Store(id='store-prev-datarate', data=2048),
                   dcc.Store(id='store-prev-channels', data=64),
                   dcc.Store(id='store-prev-subbands', data=8),
                   dcc.Store(id='store-obs-params', data=None),
                   dcc.Store(id='store-uv-data', data=None),
                   dcc.Store(id='store-elev-data', data=None),
                   dcc.Store(id='store-worldmap-data', data=None),
                   layout.top_banner(app),
                   # inputs.modal_welcome(),
                   html.Div(id='main-window', className='container-fluid d-flex row p-0 m-0',
                            children=[html.Div(id='right-column', className='col-12 col-lg-6 m-0 p-0',
                                               children=layout.inputs_column(app)),
                                      html.Div(id='left-column', className='col-12 col-lg-6 m-0 p-0',
                                               children=[layout.compute_buttons(app),
                                                         layout.outputs_column(app)])]),
                   html.Div(html.A(html.I(className="fa-solid fa-circle-info",
                    style={"font-size": "4rem"}),
                                   id="more-info-button",
                                   className="btn-floating-info btn-lg rounded-circle")),
                   dbc.Tooltip("Opens more information", target='more-info-button'),
                   dbc.Offcanvas(children=inputs.modal_general_info(), id='more-info-modal',
                                 is_open=False, className='shadow-lg blur', placement='end'),
                   html.Div(id='bottom-banner', children=[html.Br(), html.Br(), html.Br()])]))

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

@app.callback(
    Output('url', 'search', allow_duplicate=True),
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
    return search, search

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

@app.callback(
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


def main(debug: bool = False, host: str = '127.0.0.1', port: int = 8050):
    """Start the EVN Observation Planner web server.

    Parameters
    ----------
    debug : bool
        Enable Dash debug mode.
    host : str
        Host address to bind to.
    port : int
        Port number to listen on.
    """
    return app.run(debug=debug, host=host, port=port)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="EVN Observation Planner (GUI)", prog="planobs-server")
    parser.add_argument('-d', '--debug', action='store_true', default=False, help="Enable debug mode")
    parser.add_argument('--host', type=str, default='127.0.0.1', help="Host address (default: 127.0.0.1)")
    parser.add_argument('--port', type=int, default=8050, help="Port number (default: 8050)")
    args = parser.parse_args()
    main(debug=args.debug, host=args.host, port=args.port)
