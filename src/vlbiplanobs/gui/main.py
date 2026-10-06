"""Real-time version of main.py - updates outputs as inputs change (no compute button needed)."""
from __future__ import annotations
import os
import argparse
import hashlib
import threading
import time
from collections import OrderedDict
from importlib.util import find_spec
from typing import Optional
from loguru import logger
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime as dt
import flask
from plotly.io.json import to_json_plotly
from dash import Dash, html, dcc, Output, Input, State, MATCH, ALL, no_update
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from astropy.utils.iers import conf as iers_conf
from astropy import units as u

# The IERS configuration must be set before importing any vlbiplanobs module: observation.py loads the
# Earth-rotation (IERS) table when imported. If these options changed afterwards, astropy would try to
# download the table during the import, and then read the bundled table again on the first user request
# (once per computing thread: that alone used to take several seconds on the first request of each worker).
iers_conf.auto_download = False
iers_conf.auto_max_age = None

from vlbiplanobs import sources  # noqa: E402
from vlbiplanobs import observation  # noqa: E402
from vlbiplanobs import cli  # noqa: E402
from vlbiplanobs.gui import inputs, outputs, validation  # noqa: E402
from vlbiplanobs.gui.callbacks import *  # noqa: E402,F401,F403
from vlbiplanobs.gui import layout  # noqa: E402


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
# All styling assets are served locally (no third-party CDNs) so the app makes zero
# external font/CSS requests and therefore embeds no trackers (a GDPR concern for the
# previously-used Google Fonts / Font Awesome / jsDelivr CDNs). Vendored copies live
# under assets/vendor/ and the fonts under assets/css/local-fonts.css:
#   - vendor/flatly/bootstrap.min.css   : Bootswatch Flatly theme (was dbc.themes.FLATLY CDN).
#   - vendor/fontawesome/css/all.min.css : Font Awesome 6.7.2 Free (was dbc.icons.FONT_AWESOME
#                                          + the cdnjs FA 5.15.4 sheet; covers both the
#                                          "fa-solid ..." and legacy "fa fa-..." class syntaxes).
# Google Fonts @import lines were stripped from the vendored CSS; Inter/Lato now load from
# local files via css/local-fonts.css. dbc.icons.BOOTSTRAP was dropped (no "bi-" classes are
# used) and dmc.styles.DATES resolves to an empty string in this dmc version.
# Order matters: the Bootstrap theme must load before the component/soft-ui overrides.
external_stylesheets: list = ['/assets/vendor/flatly/bootstrap.min.css',
                              '/assets/vendor/fontawesome/css/all.min.css']
external_scripts: list = []

# Dash auto-loads every file under assets/. Exclude files that would otherwise throw
# console errors or double-load on every page:
#   - soft-ui-dashboard.min.js: duplicate of soft-ui-dashboard.js -> "duplicate variable 'className'".
#   - Chart.extension.js + chartjs.min.js: Chart.js plugin/lib; app charts with Plotly, not Chart.js.
#   - bootstrap-notify.js: jQuery notification plugin; jQuery is never loaded and never used.
#   - font-awesome.min.css: old Font Awesome 4.7.0 sheet whose webfonts are absent (would 404) and
#     whose ".fa" rule would clash with the vendored FA6; superseded by vendor/fontawesome.
#   - bootstrap.min.css / all.min.css: the vendored Flatly + FA6 sheets, loaded (in order) via
#     external_stylesheets above, so they must NOT also be auto-loaded from assets/ (double load).
assets_ignore = (r'soft-ui-dashboard\.min\.js|Chart\.extension\.js|chartjs\.min\.js|'
                 r'bootstrap-notify\.js|font-awesome\.min\.css|bootstrap\.min\.css|all\.min\.css')


class PlanobsDash(Dash):
    """Dash app that serializes its (static) layout only once.

    Dash builds the JSON of the layout again on every page load. Here the layout never changes and it is
    large (about 0.5 MB, mostly the cards of all the antennas), which took ~0.25 s per page load. The
    response also carries an ETag, so a browser that already has the layout only gets back an empty 304.
    """
    _layout_body: Optional[bytes] = None
    _layout_etag: str = ''

    def serve_layout(self):
        if self._layout_body is None:
            self._layout_body = super().serve_layout().get_data()
            self._layout_etag = hashlib.sha1(self._layout_body).hexdigest()

        # A substring match, as flask-compress appends ':<algorithm>' to the ETag of the compressed responses
        if self._layout_etag in flask.request.headers.get('If-None-Match', ''):
            response = flask.Response(status=304)
        else:
            response = flask.Response(self._layout_body, mimetype='application/json')

        response.set_etag(self._layout_etag)
        response.cache_control.no_cache = True
        return response


# Compress the responses (layout, callback outputs, JS/CSS bundles are 5-25 times smaller) when
# flask-compress is available (installed with the 'dash[compress]' dependency).
app = PlanobsDash(__name__, title='EVN Observation Planner', external_scripts=external_scripts,
                  external_stylesheets=external_stylesheets,
                  assets_folder=current_directory+'/assets/', assets_ignore=assets_ignore,
                  serve_locally=True,
                  eager_loading=False,
                  compress=find_spec('flask_compress') is not None,
                  suppress_callback_exceptions=True,
                  prevent_initial_callbacks=False)  # Allow initial callbacks for real-time updates


@app.server.after_request
def _cache_static_assets(response: flask.Response) -> flask.Response:
    """Let the browsers cache the files under assets/ instead of asking for every one on each page load.

    Flask serves them with 'Cache-Control: no-cache', so each of them costs a request on every visit.
    The CSS/JS files that Dash loads carry their modification time in the URL (?m=...), so they can be
    cached "forever". Everything else (images, fonts, vendored CSS) is cached for a day.
    """
    if response.status_code in (200, 304) and flask.request.path.startswith(_ASSETS_URL_PATH):
        response.cache_control.no_cache = None
        response.cache_control.public = True
        if 'm' in flask.request.args:
            response.cache_control.max_age = 31536000
            response.cache_control.immutable = True
        else:
            response.cache_control.max_age = 86400

    return response


_ASSETS_URL_PATH: str = f"{app.config.routes_pathname_prefix}{app.config.assets_url_path.strip('/')}/"


# --------------------------------------------------------------------------------------
# Per-target PDF download (one button per output tab + one for the no-target panel).
# --------------------------------------------------------------------------------------
def _params_to_obs(obs_params: dict, target_spec: Optional[str] = None) -> Optional[cli.VLBIObs]:
    """Rebuild a VLBIObs from validated observation parameters.

    The compute callback stores all the inputs it used in store-obs-params so the PDF
    download callback can reconstruct the same observation without reading the GUI state
    again. When target_spec is given, only that single target is included in the
    rebuilt observation.

    Parameters
    ----------
    obs_params : dict
        Observation parameters already normalized by validation.normalize_obs_params
        (store-obs-params is client-controlled and must never be used unvalidated).
    target_spec : str or None, optional
        Target specification to include. If None, all targets are included.

    Returns
    -------
    VLBIObs or None
        Reconstructed observation object, or None if params is None.
    """
    if obs_params is None:
        return None
    targets = [target_spec] if target_spec is not None else obs_params.get('targets')
    duration = obs_params['duration'] * u.h if obs_params['duration'] is not None else None
    return cli.main(band=obs_params['band'], stations=obs_params['stations'], targets=targets, duration=duration,
                    ontarget=obs_params['ontarget'], start_time=obs_params['start_time'],
                    datarate=obs_params['datarate'] * u.Mbit / u.s, subbands=obs_params['subbands'],
                    channels=obs_params['channels'], polarizations=obs_params['polarizations'],
                    inttime=obs_params['inttime'] * u.s)


@app.callback(Output('download-data', 'data'),
              Input('button-download', 'n_clicks'),
              State('store-obs-params', 'data'),
              prevent_initial_call=True)
def download_pdf(n_clicks: int, obs_params: dict):
    """Generate one PDF summary with a separate section for every target source.

    Parameters
    ----------
    n_clicks : int
        Number of clicks on the download button.
    obs_params : dict
        Serialized observation parameters from store-obs-params.

    Returns
    -------
    dict
        Download data for dcc.send_file.

    Raises
    ------
    PreventUpdate
        If no clicks or invalid observation parameters are available.
    """
    if not n_clicks or obs_params is None:
        raise PreventUpdate

    try:
        obs_params = validation.normalize_obs_params(obs_params)
    except validation.InvalidObsParams as e:
        logger.warning(f"PDF download rejected: invalid store-obs-params ({e}).")
        raise PreventUpdate

    try:
        target_specs = obs_params.get('targets') or [None]
        observations = [_params_to_obs(obs_params, target_spec=target_spec) for target_spec in target_specs]
        observations = [observation for observation in observations if observation is not None]
        if not observations:
            raise PreventUpdate
        logger.info(f"PDF generation started for {len(observations)} target page(s).")
        try:
            tmpfile = outputs.summary_pdf_for_sources(observations, show_figure=True)
        except Exception as fig_error:
            logger.warning(f"Could not include figures in PDF: {fig_error}")
            tmpfile = outputs.summary_pdf_for_sources(observations, show_figure=False)

        logger.info(f"PDF created at {tmpfile}.")
        download_data = dcc.send_file(tmpfile, filename='planobs_summary.pdf')
        try:
            os.remove(tmpfile)
        except OSError:
            logger.warning(f"Could not remove temporary PDF file {tmpfile}.")
        return download_data
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

    Parameters
    ----------
    target_spec : str or None
        Target specification, or None for no target.
    shared_kwargs : dict
        Shared keyword arguments for cli.main.

    Returns
    -------
    tuple[VLBIObs or None, str or None]
        (obs, error_message). error_message is None on success.
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


# Content of the output tab already built for a (target, observation setup). The real-time callback
# builds all tabs again on every change of the inputs, so when the user adds or removes a target, or gets
# back to a previous setup, the tabs of the unchanged targets are reused instead of computed again.
_TAB_CACHE: OrderedDict = OrderedDict()
_TAB_CACHE_LOCK = threading.Lock()
_TAB_CACHE_SIZE = 48
_TAB_CACHE_MAX_AGE = 3600.0  # seconds


def _tab_cache_key(target_spec: str, shared_kwargs: dict) -> tuple:
    """Return a hashable key identifying the output tab of a target under the given observation setup.

    Parameters
    ----------
    target_spec : str
        Target specification.
    shared_kwargs : dict
        Shared keyword arguments for cli.main.

    Returns
    -------
    tuple
        Hashable key. It includes the current date as the outputs without a defined epoch refer to the
        current year.
    """
    return (target_spec, dt.now().date().isoformat()) + tuple(
        (name, tuple(value) if isinstance(value, list) else str(value))
        for name, value in sorted(shared_kwargs.items()))


def _target_tab_content(target_spec: str, shared_kwargs: dict, executor: ThreadPoolExecutor):
    """Return a function providing the content of the output tab of a target, computed in the executor.

    The computation is skipped when the same tab has been built recently (see _TAB_CACHE). Tabs showing
    an error are never cached, as the error can be transient (e.g. resolving the name of the source).

    Parameters
    ----------
    target_spec : str
        Target specification.
    shared_kwargs : dict
        Shared keyword arguments for cli.main.
    executor : ThreadPoolExecutor
        Executor that runs the computation of the observation.

    Returns
    -------
    Callable
        Function without arguments that returns the content of the tab.
    """
    key = _tab_cache_key(target_spec, shared_kwargs)
    with _TAB_CACHE_LOCK:
        cached = _TAB_CACHE.get(key)
        if cached is not None and time.monotonic() - cached[0] < _TAB_CACHE_MAX_AGE:
            _TAB_CACHE.move_to_end(key)
            return lambda: cached[1]

    future = executor.submit(_compute_one_target, target_spec, shared_kwargs)

    def build():
        obs_for_target, err = future.result()
        content = outputs.build_target_tab_content(obs_for_target, target_spec, error=err)
        if obs_for_target is not None and err is None:
            with _TAB_CACHE_LOCK:
                _TAB_CACHE[key] = (time.monotonic(), content)
                _TAB_CACHE.move_to_end(key)
                while len(_TAB_CACHE) > _TAB_CACHE_SIZE:
                    _TAB_CACHE.popitem(last=False)

        return content

    return build


def _invalid_outputs(hidden_outputs: tuple, message: str) -> tuple:
    """Return the compute-callback outputs used when a user input fails server-side validation.

    Parameters
    ----------
    hidden_outputs : tuple
        The default "nothing to show" outputs of compute_observation_realtime.
    message : str
        Human-readable reason shown to the user in the user-message area.

    Returns
    -------
    tuple
        hidden_outputs with the user message replaced by an error card.
    """
    return (outputs.error_card("Invalid observation parameters", message),) + tuple(hidden_outputs[1:])


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
              [State({'type': 'network-switch', 'index': ALL}, 'value')],
              suppress_callback_exceptions=True)
def compute_observation_realtime(band: int, target_specs: Optional[list[str]], onsourcetime: int,
                                 defined_epoch: bool, startdate: str, starttime: str, duration: int | float,
                                 datarate: int, subbands: int, channels: int, pols: int, inttime: int, e_evn: bool,
                                 selected_antennas: list[str], network_switches: Optional[list[bool]] = None):
    """Real-time computation: builds the outputs container with one tab per target.

    The container shows:
    - nothing while the inputs are not enough to run any computation;
    - a single panel with duration-only outputs when no target source is specified;
    - a dbc.Tabs (one tab per target) otherwise.

    Parameters
    ----------
    band : int
        Selected band index.
    target_specs : list[str] or None
        List of target specifications.
    onsourcetime : int
        Percentage of on-source time.
    defined_epoch : bool
        Whether an epoch is specified.
    startdate : str
        Start date string.
    starttime : str
        Start time string.
    duration : int or float
        Observation duration in hours.
    datarate : int
        Data rate in Mbit/s.
    subbands : int
        Number of subbands.
    channels : int
        Number of channels.
    pols : int
        Number of polarizations.
    inttime : int
        Integration time in seconds.
    e_evn : bool
        Whether e-EVN mode is enabled.
    selected_antennas : list[str]
        List of selected antenna codenames.
    network_switches : list[bool] or None
        On/off state of each network switch, ordered as observation._NETWORKS.

    Returns
    -------
    tuple
        (user_message, loading_div, outputs_container, prev_datarate, prev_channels,
         prev_subbands, obs_params).
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

    try:
        band_name = validation.band_name_from_index(band)
    except validation.InvalidObsParams as e:
        logger.warning(f"Real-time update rejected: {e}")
        return _invalid_outputs(hidden_outputs, str(e))
    if band_name is None or not selected_antennas:
        return hidden_outputs

    selected_antennas = [ant for ant in validation.clean_station_codenames(selected_antennas)
                         if observation._STATIONS[ant].has_band(band_name)
                         and (not e_evn or observation._STATIONS[ant].real_time)]
    if not selected_antennas:
        return hidden_outputs

    target_specs = validation.clean_target_specs(target_specs)
    has_targets = bool(target_specs)
    try:
        duration = validation.validate_duration(duration)
        ontarget = validation.ontarget_from_percent(onsourcetime)
        datarate = validation.validate_datarate(datarate)
        subbands = validation.validate_subbands(subbands)
        channels = validation.validate_channels(channels)
        pols = validation.validate_polarizations(pols)
        inttime = validation.validate_inttime(inttime)
    except validation.InvalidObsParams as e:
        logger.warning(f"Real-time update rejected: {e}")
        return _invalid_outputs(hidden_outputs, str(e))
    has_duration = duration is not None

    if defined_epoch:
        epoch_complete = startdate is not None and starttime is not None and has_duration
        if not epoch_complete and has_targets:
            defined_epoch = False
        elif not epoch_complete and not has_targets and not has_duration:
            return hidden_outputs

    if not has_targets and not has_duration:
        return hidden_outputs

    start_time = None
    if defined_epoch and startdate and starttime:
        try:
            start_time = validation.parse_start_time(startdate, starttime)
        except validation.InvalidObsParams as e:
            logger.warning(f"Real-time update rejected: {e}")
            return _invalid_outputs(hidden_outputs, str(e))

    t0 = dt.now()
    shared_kwargs = dict(band=band_name, stations=sorted(selected_antennas),
                         duration=duration * u.h if has_duration else None, ontarget=ontarget, start_time=start_time,
                         datarate=datarate * u.Mbit / u.s, subbands=subbands, channels=channels,
                         polarizations=pols, inttime=inttime * u.s)

    selected_networks = [name for name, on in zip(observation._NETWORKS,
                                                   network_switches or []) if on]
    logger.info(f"Real-time update: band={band_name}, "
                f"antennas={','.join(sorted(selected_antennas))}, "
                f"networks={','.join(selected_networks) if selected_networks else 'none'}, "
                f"targets={target_specs if has_targets else 'none'}, duration={duration}")

    if has_targets:
        # Compute one observation per target in parallel (only the ones not built recently).
        with ThreadPoolExecutor(max_workers=max(1, min(len(target_specs), 8))) as executor:
            contents = {t: _target_tab_content(t, shared_kwargs, executor) for t in target_specs}

        tabs = []
        target_count = 0
        for spec, build_content in contents.items():
            # If the spec looks like coordinates (contains ':' or h/m/s pattern),
            # label as "Target N". Otherwise use the spec as-is (it's a source name).
            is_coords = ':' in spec or all(c in spec for c in ('h', 'm', 's'))
            if is_coords:
                target_count += 1
                label = f"Target {target_count}"
            else:
                label = spec
            tabs.append(dbc.Tab(build_content(), label=label, tab_id=f"tab-{spec}"))

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

    obs_params = {'band': band_name, 'stations': sorted(selected_antennas), 'targets': target_specs or None,
                  'duration': duration, 'ontarget': ontarget,
                  'startdate': startdate if start_time is not None else None,
                  'starttime': starttime if start_time is not None else None,
                  'datarate': datarate, 'subbands': subbands, 'channels': channels, 'polarizations': pols,
                  'inttime': inttime}

    return (
        html.Div(),                # user-message (empty: details are inside the tabs/panel)
        html.Div(),                # loading-div
        container_children,        # outputs-container
        datarate,                  # store-prev-datarate
        channels,                  # store-prev-channels
        subbands,                  # store-prev-subbands
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
                   dcc.Store(id='suppress-network-antenna-update', data=False),
                   # Hidden compute button (needed for callback compatibility but not shown)
                   html.Div(html.Button(id='compute-observation', style={'display': 'none'})),
                   layout.top_banner(app),
                   html.Div(id='main-window', className='container-fluid d-flex row p-0 m-0',
                            children=[html.Div(id='right-column', className='col-12 col-sm-6 m-0 p-0',
                                               children=layout.inputs_column(app)),
                                      html.Div(id='left-column', className='col-12 col-sm-6 m-0 p-0',
                                               children=[layout.compute_buttons_realtime(app),
                                                         layout.export_button_div(),
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


def warm_up() -> None:
    """Run a small observation through the whole compute and render path.

    The first observation computed by a process is much slower than the following ones (astropy, astroplan
    and plotly set up their internal tables, frame transformation graphs, templates, etc. on first use).
    Running one when the app is loaded moves that cost from the first user request to the server startup.
    With gunicorn's --preload this happens only once, in the master process, before forking the workers.

    It is skipped when the environment variable PLANOBS_NO_WARMUP is set. A failure here is only logged:
    it must never prevent the server from starting.
    """
    if os.environ.get('PLANOBS_NO_WARMUP'):
        return

    t0 = time.perf_counter()
    try:
        band = '18cm'
        stations = [ant.codename for ant in observation._STATIONS if ant.has_band(band)][:6]
        shared_kwargs = dict(band=band, stations=sorted(stations), ontarget=0.7, datarate=1024 * u.Mbit / u.s,
                             subbands=8, channels=64, polarizations=4, inttime=2 * u.s)
        # Coordinates instead of a source name, so nothing needs to be resolved online. Both the default
        # (no epoch) and the fixed-epoch modes, as they follow different code paths.
        for target, duration, start_time in (('12h29m06.7s +02d03m08.6s', None, None),
                                             ('12h29m06.7s +02d03m08.6s', 4 * u.h,
                                              validation.parse_start_time(dt.now().date().isoformat(), '12:00')),
                                             (None, 4 * u.h, None)):
            obs, err = _compute_one_target(target, dict(shared_kwargs, duration=duration, start_time=start_time))
            if obs is None:
                logger.warning(f"Warm-up observation could not be computed: {err}")
            elif target is not None:
                to_json_plotly(outputs.build_target_tab_content(obs, target))
            else:
                outputs.build_no_target_panel(obs)

        # Serializes (and caches) the layout, so the workers forked afterwards share it
        with app.server.test_request_context('/_dash-layout'):
            app.serve_layout()

        logger.info(f"Warm-up completed in {time.perf_counter() - t0:.2f}s")
    except Exception as e:
        logger.warning(f"Warm-up failed ({type(e).__name__}: {e}); the first request will be slower.")


warm_up()


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
