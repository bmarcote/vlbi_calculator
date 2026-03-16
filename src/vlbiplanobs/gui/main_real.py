"""Real-time version of main.py - updates outputs as inputs change (no compute button needed)."""
import os
import argparse
from loguru import logger
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime as dt
from dash import Dash, html, dcc, Output, Input, State, no_update
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

app = Dash(__name__, title='EVN Observation Planner (Real-time)', external_scripts=external_scripts,
           external_stylesheets=[dbc.themes.FLATLY, dbc.icons.BOOTSTRAP,
                                 dbc.icons.FONT_AWESOME, dmc.styles.DATES] + external_stylesheets,
           assets_folder=current_directory+'/assets/', eager_loading=False,
           prevent_initial_callbacks=False)  # Allow initial callbacks for real-time updates


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


# Real-time callback - triggers on ANY input change
# Plot figures are rendered progressively via store-triggered callbacks in callbacks.py
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
              [Input('band-slider', 'value'),
               Input('switch-specify-source', 'value'),
               Input('source-input', 'value'),
               Input('onsourcetime', 'value'),
               Input('switch-specify-epoch', 'value'),
               Input('startdate', 'value'),
               Input('starttime', 'value'),
               Input('duration', 'value'),
               Input('datarate', 'value'),
               Input('subbands', 'value'),
               Input('channels', 'value'),
               Input('pols', 'value'),
               Input('inttime', 'value'),
               Input('switch-specify-e-evn', 'value'),
               Input('switches-antennas', 'value')],
              suppress_callback_exceptions=True)
def compute_observation_realtime(band: int, defined_source: bool, source: str, onsourcetime: int,
                                  defined_epoch: bool, startdate: str, starttime: str, 
                                  duration: int | float, datarate: int, subbands: int,
                                  channels: int, pols: int, inttime: int, e_evn: bool, 
                                  selected_antennas: list[str]):
    """Real-time computation - returns text outputs immediately, plots render via stores.

    Shows whatever outputs are possible with the current inputs:
    - No band/antennas: Show nothing
    - Band + antennas + (source OR duration): Show full observation
    """
    # Default hidden state for all outputs (27 values)
    hidden_outputs = (
        html.Div(),  # user-message
        html.Div(),  # loading-div
        True,        # download-summary-div hidden
        no_update,   # download-link
        html.Div(),  # card-rms
        html.Div(),  # sensitivity-baseline-modal
        html.Div(),  # card-resolution
        html.Div(),  # out-sun
        html.Div(),  # out-phaseref
        html.Div(),  # out-ant
        True,        # out-elevations hidden
        html.Div(),  # out-elevations-info
        True,        # out-uv-coverage hidden
        html.Div(),  # out-uv-coverage-info
        [],          # select-antenna-uv-plot options
        html.Div(),  # div-card-fov
        html.Div(),  # div-card-vel
        True,        # out-worldmap hidden
        html.Div(),  # card-datasize
        html.Div(),  # div-card-time
        no_update,   # store-prev-datarate
        no_update,   # store-prev-channels
        no_update,   # store-prev-subbands
        None,        # store-obs-params
        None,        # store-uv-data
        None,        # store-elev-data
        None,        # store-worldmap-data
    )
    
    if band == 0 or band is None:
        return hidden_outputs
    
    if not selected_antennas:
        return hidden_outputs
    
    selected_antennas = [ant for ant in selected_antennas
                         if observation._STATIONS[ant].has_band(inputs.band_from_index(band))
                         and (not e_evn or observation._STATIONS[ant].real_time)]
    
    if not selected_antennas:
        return hidden_outputs
    
    has_source = defined_source and source and source.strip() != ''
    has_duration = duration is not None and duration >= 0.1
    
    if defined_epoch:
        epoch_complete = startdate is not None and starttime is not None and has_duration
        if not epoch_complete and has_source:
            defined_epoch = False
        elif not epoch_complete and not has_source:
            if not has_duration:
                return hidden_outputs
    
    if not has_source and not has_duration:
        return hidden_outputs
    
    t0 = dt.now()
    try:
        logger.info(f"Real-time update: band={inputs.band_from_index(band)}, "
                    f"antennas={len(selected_antennas)}, source={source}, duration={duration}")
        
        obs = cli.main(
            band=inputs.band_from_index(band),
            stations=sorted(selected_antennas),
            targets=[source] if has_source else None,
            duration=duration * u.h if has_duration else None,
            ontarget=onsourcetime / 100 if onsourcetime else 0.7,
            start_time=Time(dt.strptime(f"{startdate} {starttime}", '%Y-%m-%d %H:%M'),
                           format='datetime', scale='utc') if defined_epoch and startdate and starttime else None,
            datarate=(datarate if not isinstance(datarate, str) else int(datarate or 2048)) * u.Mbit / u.s,
            subbands=subbands or 8,
            channels=channels or 64,
            polarizations=pols or 2,
            inttime=(inttime or 2) * u.s
        )
        
        if obs is None:
            return hidden_outputs
        
        with ThreadPoolExecutor() as executor:
            f_observable = executor.submit(obs.is_observable)
            f_sun_constraint = executor.submit(obs.sun_constraint)
            f_uv_data = executor.submit(obs.get_uv_data)
            f_baseline_sens = executor.submit(obs.baseline_sensitivity)
            _ = f_observable.result()
            _ = f_sun_constraint.result()
            _ = f_uv_data.result()
            _ = f_baseline_sens.result()

        with ThreadPoolExecutor() as executor:
            f_always_obs = executor.submit(obs.is_always_observable)
            f_sun_epochs = executor.submit(obs.sun_limiting_epochs)
            f_thermal = executor.submit(obs.thermal_noise)
            f_uv_values = executor.submit(obs.get_uv_values)
            _ = f_always_obs.result()
            _ = f_sun_epochs.result()
            _ = f_thermal.result()
            _ = f_uv_values.result()

        _ = obs.synthesized_beam()
        
    except (ValueError, sources.SourceNotVisible, IndexError) as e:
        logger.debug(f"Cannot compute observation yet: {e}")
        return hidden_outputs
    except Exception as e:
        logger.debug(f"Error during real-time update: {e}")
        return hidden_outputs
    
    try:
        if not all(obs.is_observable_by_network(min_stations=1).values()):
            return hidden_outputs
        
        futures = {}
        can_show_source_plots = obs.sourcenames and has_source
        
        with ThreadPoolExecutor() as executor:
            # Text outputs (fast)
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
            futures['worldmap_data'] = executor.submit(plots.serialize_worldmap_data, obs)
            
            if can_show_source_plots:
                futures['out_sun'] = executor.submit(outputs.sun_warning, obs)
                futures['elev_data'] = executor.submit(plots.serialize_elevation_data, obs)
                futures['uv_data'] = executor.submit(plots.serialize_uv_data, obs)
                futures['out_elev_info'] = executor.submit(outputs.print_observability_ranges, obs)
                futures['out_uv_info'] = executor.submit(outputs.print_baseline_lengths, obs)
                futures['out_ant_options'] = executor.submit(outputs.put_antenna_options, obs)

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

            if can_show_source_plots:
                out_sun = futures['out_sun'].result()
                out_elev = [False, futures['out_elev_info'].result()]
                out_uv = [False, futures['out_uv_info'].result(), futures['out_ant_options'].result()]
                uv_data = futures['uv_data'].result()
                elev_data = futures['elev_data'].result()
            else:
                out_sun = html.Div()
                out_elev = [True, html.Div()]
                out_uv = [True, html.Div(), []]
                uv_data = None
                elev_data = None
                
    except Exception as e:
        logger.debug(f"Error generating outputs: {e}")
        return hidden_outputs

    elapsed = (dt.now() - t0).total_seconds()
    logger.info(f"Real-time update completed in {elapsed:.2f}s")

    obs_params = {
        'band': inputs.band_from_index(band),
        'stations': sorted(selected_antennas),
        'targets': [source] if has_source else None,
        'duration': duration,
        'ontarget': onsourcetime / 100 if onsourcetime else 0.7,
        'startdate': startdate if defined_epoch else None,
        'starttime': starttime if defined_epoch else None,
        'datarate': datarate if not isinstance(datarate, str) else int(datarate or 2048),
        'subbands': subbands or 8,
        'channels': channels or 64,
        'polarizations': pols or 2,
        'inttime': inttime or 2
    }

    return (
        html.Div(),        # user-message
        html.Div(),        # loading-div
        False,             # download-summary-div
        None,              # download-link
        out_rms,           # card-rms
        out_baseline_sens, # sensitivity-baseline-modal
        out_res,           # card-resolution
        out_sun,           # out-sun
        out_phaseref,      # out-phaseref
        out_ant,           # out-ant
        *out_elev,         # out-elevations (hidden, info)
        *out_uv,           # out-uv-coverage (hidden, info, options)
        out_fov,           # div-card-fov
        out_freq,          # div-card-vel
        False,             # out-worldmap
        out_datasize,      # card-datasize
        out_obstime,       # div-card-time
        datarate or 2048,  # store-prev-datarate
        channels or 64,    # store-prev-channels
        subbands or 8,     # store-prev-subbands
        obs_params,        # store-obs-params
        uv_data,           # store-uv-data
        elev_data,         # store-elev-data
        worldmap_data,     # store-worldmap-data
    )


server = app.server
app.index_string = app.index_string.replace('<body>', '<body class="g-sidenav-show bg-gray-100">')

# Layout without the compute button prominently featured
app.layout = dmc.MantineProvider(dbc.Container(fluid=True, className='bg-gray-100 row m-0 p-4', children=[
                   dcc.Store(id='store-prev-datarate', data=2048),
                   dcc.Store(id='store-prev-channels', data=64),
                   dcc.Store(id='store-prev-subbands', data=8),
                   dcc.Store(id='store-obs-params', data=None),
                   dcc.Store(id='store-uv-data', data=None),
                   dcc.Store(id='store-elev-data', data=None),
                   dcc.Store(id='store-worldmap-data', data=None),
                   # Hidden compute button (needed for callback compatibility but not shown)
                   html.Div(html.Button(id='compute-observation', style={'display': 'none'})),
                   layout.top_banner(app),
                   html.Div(id='main-window', className='container-fluid d-flex row p-0 m-0',
                            children=[html.Div(id='right-column', className='col-12 col-sm-6 m-0 p-0',
                                               children=layout.inputs_column(app)),
                                      html.Div(id='left-column', className='col-12 col-sm-6 m-0 p-0',
                                               children=[layout.compute_buttons_realtime(app),
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
    description = "EVN Observation Planner (Real-time GUI)"
    parser = argparse.ArgumentParser(description=description, prog="planobs-server-realtime", usage=usage)
    parser.add_argument('-d', '--debug', action='store_true', default=False, help="Enable debug mode")
    args = parser.parse_args()
    return app.run(debug=args.debug)


if __name__ == '__main__':
    main(debug=False)
