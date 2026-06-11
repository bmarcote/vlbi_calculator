"""Real-time version of main.py - updates outputs as inputs change (no compute button needed)."""
from __future__ import annotations
import os
import argparse
from typing import Optional
from loguru import logger
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime as dt
from dash import Dash, html, dcc, Output, Input, State, MATCH, no_update
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from astropy.utils.iers import conf as iers_conf
from astropy import units as u
from astropy.time import Time
from vlbiplanobs import sources
from vlbiplanobs import observation
from vlbiplanobs import cli
from vlbiplanobs.gui import inputs, outputs
from vlbiplanobs.gui.callbacks import *  # noqa: F401,F403
from vlbiplanobs.gui import layout


iers_conf.auto_download = False
iers_conf.auto_max_age = None

def setup_file_logging(logfilename: Optional[str] = None) -> int:
    """Enable loguru file logging for planobs and return the sink handler id.

    No log file is created unless this function is called explicitly (e.g. via
    the CLI ``--logging`` flag). This avoids creating a log file by default.

    Parameters
    ----------
    logfilename : Optional[str]
        Path of the log file. When None, '/var/log/planobs.log' is used if
        writable, otherwise '~/log-planobs.log'.

    Returns
    -------
    int
        The loguru handler id of the added file sink.
    """
    if logfilename is None:
        if os.access("/var/log/planobs.log", os.W_OK):
            logfilename = "/var/log/planobs.log"
        else:
            logfilename = os.path.expanduser("~/log-planobs.log")

    return logger.add(logfilename, backtrace=True, diagnose=True,
                      format="{time:YYYY-MM-DD HH:mm} |  {level} {message}")

current_directory = os.path.dirname(os.path.realpath(__file__))
external_stylesheets: list = ['https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css']
external_scripts: list = []

app = Dash(__name__, title='EVN Observation Planner', external_scripts=external_scripts,
           external_stylesheets=[dbc.themes.FLATLY, dbc.icons.BOOTSTRAP,
                                 dbc.icons.FONT_AWESOME, dmc.styles.DATES] + external_stylesheets,
           assets_folder=current_directory+'/assets/', eager_loading=False,
           suppress_callback_exceptions=True,
           prevent_initial_callbacks=False)  # Allow initial callbacks for real-time updates


# --------------------------------------------------------------------------------------
# Per-target PDF download (one button per output tab + one for the no-target panel).
# --------------------------------------------------------------------------------------
def _params_to_obs(obs_params: dict, target_spec: Optional[str] = None) -> Optional[cli.VLBIObs]:
    """Rebuild a VLBIObs from serialised observation parameters.

    The compute callback stores all the inputs it used in `store-obs-params` so the PDF
    download callback can reconstruct the same observation without reading the GUI state
    again. When ``target_spec`` is given, only that single target is included in the
    rebuilt observation; otherwise the rebuild matches what the user saw on the page.
    """
    if obs_params is None:
        return None
    targets = [target_spec] if target_spec is not None else obs_params.get('targets')
    return cli.main(
        band=obs_params['band'],
        stations=obs_params['stations'],
        targets=targets,
        duration=obs_params['duration'] * u.h if obs_params['duration'] is not None else None,
        ontarget=obs_params['ontarget'],
        start_time=Time(dt.strptime(f"{obs_params['startdate']} {obs_params['starttime']}",
                                    '%Y-%m-%d %H:%M'),
                        format='datetime', scale='utc') if obs_params.get('startdate') else None,
        datarate=obs_params['datarate'] * u.Mbit / u.s,
        subbands=obs_params['subbands'],
        channels=obs_params['channels'],
        polarizations=obs_params['polarizations'],
        inttime=obs_params['inttime'] * u.s)


@app.callback(Output({'type': 'download-pdf', 'index': MATCH}, 'data'),
              Input({'type': 'btn-pdf', 'index': MATCH}, 'n_clicks'),
              State({'type': 'btn-pdf', 'index': MATCH}, 'id'),
              State('store-obs-params', 'data'),
              prevent_initial_call=True)
def download_pdf_per_target(n_clicks, btn_id, obs_params: dict):
    """Generate and download a PDF summary for the target identified by the clicked button.

    The button id is `{'type': 'btn-pdf', 'index': <target_spec>}`. A special index
    ``'__no_target__'`` is used for the duration-only panel.
    """
    if not n_clicks or obs_params is None:
        raise PreventUpdate

    target_spec = btn_id['index']
    use_target: Optional[str] = None if target_spec == '__no_target__' else target_spec
    safe_label = 'summary' if use_target is None else target_spec
    safe_label = ''.join(ch if ch.isalnum() or ch in ('-', '_') else '_' for ch in safe_label)

    try:
        logger.info(f"PDF generation started for target='{target_spec}'.")
        obs = _params_to_obs(obs_params, target_spec=use_target)
        if obs is None:
            raise PreventUpdate
        try:
            tmpfile = outputs.summary_pdf(obs, show_figure=True)
        except Exception as fig_error:
            logger.warning(f"Could not include figure in PDF: {fig_error}")
            tmpfile = outputs.summary_pdf(obs, show_figure=False)
        logger.info(f"PDF created at {tmpfile}.")
        return dcc.send_file(tmpfile, filename=f"planobs_{safe_label}.pdf")
    except Exception as e:
        logger.exception(f"While downloading the PDF: {e}")
        raise PreventUpdate


# --------------------------------------------------------------------------------------
# Per-tab sensitivity-baseline modal toggle.
# --------------------------------------------------------------------------------------
app.clientside_callback(
    "function(n_clicks, is_open) { return n_clicks ? !is_open : dash_clientside.no_update; }",
    Output({'type': 'modal-sens', 'index': MATCH}, 'is_open'),
    Input({'type': 'btn-sens', 'index': MATCH}, 'n_clicks'),
    State({'type': 'modal-sens', 'index': MATCH}, 'is_open'),
    prevent_initial_call=True)


# --------------------------------------------------------------------------------------
# Fix for Plotly graphs in hidden tabs: trigger a window resize when user switches tabs
# so graphs re-render at correct dimensions.
# --------------------------------------------------------------------------------------
app.clientside_callback(
    """function(active_tab) {
        setTimeout(function() { window.dispatchEvent(new Event('resize')); }, 100);
        return dash_clientside.no_update;
    }""",
    Output('outputs-container', 'className'),
    Input('outputs-tabs', 'active_tab'),
    prevent_initial_call=True)


# --------------------------------------------------------------------------------------
# Real-time multi-target compute callback.
# --------------------------------------------------------------------------------------
def _compute_one_target(target_spec: Optional[str], shared_kwargs: dict) -> tuple[Optional[cli.VLBIObs], Optional[str]]:
    """Run cli.main for a single target (or no target) and warm up its caches.

    Returns ``(obs, error_message)``. ``error_message`` is None on success.
    """
    targets = [target_spec] if target_spec is not None else None
    try:
        obs = cli.main(targets=targets, **shared_kwargs)
        if obs is None:
            return None, "cli.main returned no observation."
        # Warm up the heaviest caches so all output helpers reuse them.
        obs.is_observable()
        obs.sun_constraint()
        if obs.scans:
            obs.get_uv_data()
            obs.baseline_sensitivity()
            obs.is_always_observable()
            obs.sun_limiting_epochs()
            obs.thermal_noise()
            obs.synthesized_beam()
        return obs, None
    except sources.SourceNotVisible as e:
        return None, f"Source not visible from the array: {e}"
    except ValueError as e:
        return None, str(e)
    except Exception as e:
        logger.debug(f"Error computing target '{target_spec}': {e}")
        return None, f"{type(e).__name__}: {e}"


@app.callback([Output('user-message', 'children'),
               Output('loading-div', 'children'),
               Output('outputs-container', 'children'),
               Output('store-prev-datarate', 'data'),
               Output('store-prev-channels', 'data'),
               Output('store-prev-subbands', 'data'),
               Output('store-obs-params', 'data')],
              [Input('band-slider', 'value'),
               Input('store-targets', 'data'),
               Input('onsourcetime', 'value'),
               Input('switch-specify-epoch', 'value'),
               Input('startdate', 'date'),
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
def compute_observation_realtime(band: int,
                                  target_specs: Optional[list[str]],
                                  onsourcetime: int,
                                  defined_epoch: bool, startdate: str, starttime: str,
                                  duration: int | float, datarate: int, subbands: int,
                                  channels: int, pols: int, inttime: int, e_evn: bool,
                                  selected_antennas: list[str]):
    """Real-time computation: builds the outputs container with one tab per target.

    The container shows:
    - nothing while the inputs are not enough to run any computation;
    - a single panel with duration-only outputs when no target source is specified;
    - a `dbc.Tabs` (one tab per target) otherwise.
    """
    empty_message = html.Blockquote(
        className='text-secondary text-bold ms-2 px-2',
        children=["Set your VLBI observation on the left to see the details of the "
                  "expected outcome here.",
                  html.Footer(className='text-sm pt-2',
                              children="Add target sources via the 'Source & Epoch' "
                                       "panel to compare them side-by-side.")])
    hidden_outputs = (empty_message, html.Div(), html.Div(),
                      no_update, no_update, no_update, None)

    if band == 0 or band is None or not selected_antennas:
        return hidden_outputs

    selected_antennas = [ant for ant in selected_antennas
                         if observation._STATIONS[ant].has_band(inputs.band_from_index(band))
                         and (not e_evn or observation._STATIONS[ant].real_time)]
    if not selected_antennas:
        return hidden_outputs

    target_specs = [t for t in (target_specs or []) if t and t.strip()]
    has_targets = bool(target_specs)
    has_duration = duration is not None and duration >= 0.1

    if defined_epoch:
        epoch_complete = startdate is not None and starttime is not None and has_duration
        if not epoch_complete and has_targets:
            defined_epoch = False
        elif not epoch_complete and not has_targets and not has_duration:
            return hidden_outputs

    if not has_targets and not has_duration:
        return hidden_outputs

    t0 = dt.now()
    shared_kwargs = dict(
        band=inputs.band_from_index(band),
        stations=sorted(selected_antennas),
        duration=duration * u.h if has_duration else None,
        ontarget=onsourcetime / 100 if onsourcetime else 0.7,
        start_time=Time(dt.strptime(f"{startdate[:10]} {starttime}", '%Y-%m-%d %H:%M'),
                        format='datetime', scale='utc')
        if defined_epoch and startdate and starttime else None,
        datarate=(datarate if not isinstance(datarate, str) else int(datarate or 2048)) * u.Mbit / u.s,
        subbands=subbands or 8,
        channels=channels or 64,
        polarizations=pols or 2,
        inttime=(inttime or 2) * u.s)

    logger.info(f"Real-time update: band={inputs.band_from_index(band)}, "
                f"antennas={len(selected_antennas)}, "
                f"targets={target_specs if has_targets else 'none'}, duration={duration}")

    if has_targets:
        # Compute one observation per target in parallel.
        with ThreadPoolExecutor(max_workers=max(1, min(len(target_specs), 8))) as executor:
            futures = {t: executor.submit(_compute_one_target, t, shared_kwargs)
                       for t in target_specs}
            results = {t: f.result() for t, f in futures.items()}

        tabs = []
        target_count = 0
        for spec, (obs_for_target, err) in results.items():
            # If the spec looks like coordinates (contains ':' or h/m/s pattern),
            # label as "Target N". Otherwise use the spec as-is (it's a source name).
            is_coords = ':' in spec or all(c in spec for c in ('h', 'm', 's'))
            if is_coords:
                target_count += 1
                label = f"Target {target_count}"
            else:
                label = spec
            content = outputs.build_target_tab_content(obs_for_target, spec, error=err)
            tabs.append(dbc.Tab(content, label=label, tab_id=f"tab-{spec}"))

        last_tab_id = f"tab-{target_specs[-1]}"
        container_children = html.Div(dbc.Tabs(tabs, id='outputs-tabs',
                                               active_tab=last_tab_id),
                                      key=str(len(target_specs)))
    else:
        # No target specified — duration-only panel.
        obs_only, err = _compute_one_target(None, shared_kwargs)
        if obs_only is None:
            container_children = outputs.error_card(
                "Could not compute the observation",
                err or "Unknown error.")
        else:
            container_children = outputs.build_no_target_panel(obs_only)

    elapsed = (dt.now() - t0).total_seconds()
    logger.info(f"Real-time update completed in {elapsed:.2f}s")

    obs_params = {
        'band': inputs.band_from_index(band),
        'stations': sorted(selected_antennas),
        'targets': target_specs if has_targets else None,
        'duration': duration,
        'ontarget': onsourcetime / 100 if onsourcetime else 0.7,
        'startdate': startdate if defined_epoch else None,
        'starttime': starttime if defined_epoch else None,
        'datarate': datarate if not isinstance(datarate, str) else int(datarate or 2048),
        'subbands': subbands or 8,
        'channels': channels or 64,
        'polarizations': pols or 2,
        'inttime': inttime or 2,
    }

    return (
        html.Div(),                # user-message (empty: details are inside the tabs/panel)
        html.Div(),                # loading-div
        container_children,        # outputs-container
        datarate or 2048,          # store-prev-datarate
        channels or 64,            # store-prev-channels
        subbands or 8,             # store-prev-subbands
        obs_params,                # store-obs-params
    )


server = app.server
app.index_string = app.index_string.replace('<body>', '<body class="g-sidenav-show bg-gray-100">')

# Layout without the compute button prominently featured
app.layout = dmc.MantineProvider(dbc.Container(fluid=True, className='bg-gray-100 row m-0 p-4', children=[
                   dcc.Store(id='store-prev-datarate', data=2048),
                   dcc.Store(id='store-prev-channels', data=64),
                   dcc.Store(id='store-prev-subbands', data=8),
                   dcc.Store(id='store-obs-params', data=None),
                   # Persists the user's list of target source specs (strings) across reloads.
                   dcc.Store(id='store-targets', data=[], storage_type='local'),
                   # Hidden compute button (needed for callback compatibility but not shown)
                   html.Div(html.Button(id='compute-observation', style={'display': 'none'})),
                   layout.top_banner(app),
                   html.Div(id='main-window', className='container-fluid d-flex row p-0 m-0',
                            children=[html.Div(id='right-column', className='col-12 col-sm-6 m-0 p-0',
                                               children=layout.inputs_column(app)),
                                      html.Div(id='left-column', className='col-12 col-sm-6 m-0 p-0',
                                               children=[layout.compute_buttons_realtime(app),
                                                         layout.outputs_column(app)])]),
                   # Modal allowing the user to add/remove multiple target sources.
                   inputs.target_sources_modal(),
                   html.Div(html.A(html.I(className="fa-solid fa-circle-info", style={"font-size": "4rem"}),
                                   id="more-info-button",
                                   className="btn-floating-info btn-lg rounded-circle")),
                   dbc.Tooltip("Opens more information", target='more-info-button'),
                   dbc.Offcanvas(children=inputs.modal_general_info(), id='more-info-modal',
                                 is_open=False, className='shadow-lg blur', placement='end'),
                   html.Div(id='bottom-banner', children=[html.Br(), html.Br(), html.Br()])]))


def main(debug: bool = False, host: str = '127.0.0.1', port: int = 8050,
         logging: bool | str = False):
    """Start the EVN Observation Planner web server.

    Parameters
    ----------
    debug : bool
        Enable Dash debug mode.
    host : str
        Host address to bind to.
    port : int
        Port number to listen on.
    logging : bool | str
        Enable file logging. When False (default) no log file is created.
        When True a default path is used; when a string, it is the log path.
    """
    if logging:
        setup_file_logging(None if logging is True else logging)

    return app.run(debug=debug, host=host, port=port)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="EVN Observation Planner (GUI)", prog="planobs-server")
    parser.add_argument('-d', '--debug', action='store_true', default=False, help="Enable debug mode")
    parser.add_argument('--host', type=str, default='127.0.0.1', help="Host address (default: 127.0.0.1)")
    parser.add_argument('--port', type=int, default=8050, help="Port number (default: 8050)")
    parser.add_argument('--logging', nargs='?', const=True, default=False, metavar='LOGFILE',
                        help="Enable logging to a file (disabled by default).")
    args = parser.parse_args()
    main(debug=args.debug, host=args.host, port=args.port, logging=args.logging)
