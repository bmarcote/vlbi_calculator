    #! /usr/bin/env python
# -*- coding: utf-8 -*-
"""EVN Observation Planner.

Program to compute the source elevation visibility
and expected thermal noise level for a given EVN observation.
"""

__author__ = "Benito Marcote"
__credits__ = "Benito Marcote"
__license__ = "LGPLv3+"
__date__ = "2020/10/26"
__version__ = "1.0.1"
__maintainer__ = "Benito Marcote"
__email__ = "marcote@jive.eu"
__status__ = "Production"   # Prototype, Development, Production.

import os
from os import path
from time import sleep
import itertools
from importlib import resources
import datetime
from datetime import datetime as dt
import numpy as np
import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from astropy.time import Time
from astropy import coordinates as coord
from astropy import units as u

## THIS WILL NEED TO GO AWAY IN THE NEW VERSION OF ASTROPY, WHICH IS STILL NOT
## SUPPORTED BY THE CURRENT VERSION OF ASTROPLAN
# Tweak to not let astroplan crashing...
# from astropy.utils.data import clear_download_cache
# from astropy.utils import iers
# clear_download_cache()  # to be sure it is really working
#
# iers.conf.auto_download = False
# iers.conf.iers_auto_url = None
# iers.conf.auto_max_age = None
# iers.conf.remote_timeout = 100.0
# iers.conf.download_cache_lock_attempts = 10

from astroplan import FixedTarget

current_directory = path.dirname(path.realpath(__file__))
if path.isfile(current_directory + '/.astropy/cache/download/py3/lock'):
    os.remove(current_directory + '/.astropy/cache/download/py3/lock')
#########   All the previous part will be removed with astropy 4.1+ and astroplan 0.7+

from vlbiplanobs import freqsetups as fs
from vlbiplanobs import stations
from vlbiplanobs import observation
from vlbiplanobs import graphical_elements as ge
# adding the possibility of disabled. Will be implemented in a future version of dash_bootstrap_components
from vlbiplanobs.Checkbox import Checkbox


all_antennas = stations.Stations.get_stations_from_configfile()
default_arrays = stations.Stations.get_network_names_from_configfile()


sorted_networks = {'EVN': 'EVN: European VLBI Network', 'eMERLIN': 'eMERLIN (out-stations)',
                   'VLBA': 'VLBA: Very Long Baseline Array',
                   'LBA': 'LBA: Australian Long Baseline Array',
                   'KVN': 'KVN: Korean VLBI Network',
                   'Other': 'Other antennas',
                   'Decom': 'Decommissioned antennas'}

# Safety check that all these antennas are available in the file
for a_array in default_arrays:
    for a_station in default_arrays[a_array]['default_antennas']:
        assert a_station in all_antennas.codenames

doc_files = {'About the EVN Observation Planner': 'doc-contact.md',
             'About the antennas': 'doc-antennas.md',
             'Technical background': 'doc-estimations.md'}

selected_band = '6cm'
obs = observation.Observation()

external_stylesheets = []
external_scripts = ["https://kit.fontawesome.com/69c65a0ab5.js"]


app = dash.Dash(__name__, title='EVN Observation Planner', external_scripts=external_scripts,
                assets_folder=current_directory+'/assets/')

app.config.suppress_callback_exceptions = True  # Avoids error messages for id's that haven't been loaded yet
server = app.server


def get_doc_text():
    """Reads the doc files and returns it as a Div object.
    """
    temp = []
    for i,a_topic in enumerate(doc_files):
        with resources.open_text("doc", doc_files[a_topic]) as f:
            # Some files will have references to images/files in the form '{src:filename}'
            # We parse this
            parsed_text = f.read()
            while '{src:' in parsed_text:
                i0 = parsed_text.index('{src:')
                i1 = i0 + parsed_text[i0:].index('}')
                filename = parsed_text[i0+5:i1]
                parsed_text = parsed_text.replace( parsed_text[i0:i1+1],
                                                   app.get_asset_url(filename) )

            if a_topic == 'About the antennas':
                temp += [ge.create_accordion_card(a_topic,
                   [dcc.Markdown(parsed_text), ge.antenna_cards(app, all_antennas)], id=str(i), is_open=False)]
            else:
                temp += [ge.create_accordion_card(a_topic, dcc.Markdown(parsed_text),
                         id=str(i), is_open=False)]

    return html.Div(temp, className='col-12 accordion shadow-1-strong')


@app.callback([Output(f"collapse-{i}", "is_open") for i in range(len(doc_files))],
        [Input(f"group-{i}-toggle", "n_clicks") for i in range(len(doc_files))],
        [State(f"collapse-{i}", "is_open") for i in range(len(doc_files))])
def toggle_accordion(*args):
    """Allows the expansion/collapse of an HTML accordion block.
    """
    defaults = list(args[len(doc_files):])
    ctx = dash.callback_context
    if not ctx.triggered:
        return [dash.no_update]*len(doc_files)
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]

    for i in range(len(doc_files)):
        if (button_id == f"group-{i}-toggle") and (args[i] is not None):
            defaults[i] = not defaults[i]

    return defaults


def error_text(an_error):
    """Standard error message written in a modal error window.
    It returns a str mentioning 'an_error' and the contact details to report it.
    """
    return f"An error occured.\n{an_error}.\nPlease report to marcote@jive.eu " \
            "or in https://github.com/bmarcote/vlbi_calculator."


def convert_colon_coord(colon_coord):
    """Converts some coordinates given in a str format 'HH:MM:SS DD:MM:SS' to
    'HHhMMmSSs DDdMMdSSs'.
    If ':' are not present in colon_coord, then it returns the same str.
    """
    if ':' not in colon_coord:
        return colon_coord
    for l in ('h', 'm', 'd', 'm'):
        colon_coord = colon_coord.replace(':', l, 1)

    return ' '.join([f"{s}s" for s in colon_coord.split()])


def alert_message(message, title="Warning!"):
    """Produces an alert-warning message.
    'message' can be either a string or a list with different string/dash components.
    """
    if type(message) == str:
        return [html.Br(), \
                dbc.Alert([html.H4(title, className='alert-heading'), message], \
                        color='warning', dismissable=True)]
    else:
        return [html.Br(), \
                dbc.Alert([html.H4(title, className='alert-heading'), *message], \
                        color='warning', dismissable=True)]


def update_sensitivity(obs):
    """Given the observation, it sets the text for all summary cards
    with information about the observation.
    """
    cards = []
    # The time card
    cards += ge.summary_card_times(app, obs)
    cards += ge.summary_card_frequency(app, obs)
    cards += ge.summary_card_antennas(app, obs)
    cards += ge.summary_card_beam(app, obs)
    cards += ge.summary_card_rms(app, obs)
    cards += ge.summary_card_fov(app, obs)
    return [html.Div(className='card-deck col-12 justify-content-center',
                     children=cards)]


def arrays_with_band(arrays, a_band):
    """Returns the arrays that can observe the given band with at least two antennas.
    It excludes e-EVN if it is included in arrays.

    Inputs
    - arrays : dict
        The keys are the name of the array and the values must be a dict with the codenames
        of the antennas in an array inside the value 'default_antennas'.
    - a_band : str
        The band to be observed, following the criteria in fs.bands.

    Returns
    - arrays_with_band : str
        Comma-separated list of the arrays that can observe the given band.
    """
    tmp = []   # the list of arrays that can observe the given band
    for an_array in arrays:
        if an_array != 'e-EVN':
            if a_band in arrays[an_array]['observing_bands']:
                tmp.append(an_array)

    if len(tmp) == 0:
        return 'none'
    elif len(tmp) == 2:
        return ' and '.join(tmp)
    elif len(tmp) in (1, 3):
        return ', '.join(tmp)
    else: # >= 4
        return ', '.join(tmp[:-1]) + ' and ' + tmp[-1]



@app.callback(Output('initial-pickband-label', 'children'),
              [Input('initial-band', 'value')])
def update_pickband_tooltip(a_wavelength):
    a_band = tuple(fs.bands)[a_wavelength]
    return [dbc.Card(dbc.CardBody([
                    html.H5([html.Img(height='30rem',
                                      src=app.get_asset_url(f"waves-{a_band.replace('.',  '_')}.png"),
                                      alt='Band: ', className='d-inline-block'),
                             html.Span(f"{fs.bands[a_band].split('(')[0].strip()}",
                                       style={'float': 'right'})
                             ], className="card-title"),
                    html.P([html.Span("Wavelength: ", style={'color': '#888888'}),
                        f"{fs.bands[a_band].split('(')[1].split('or')[0].strip()}.",
                        html.Br(),
                        html.Span("Frequency: ", style={'color': '#888888'}),
                        f"{fs.bands[a_band].split('(')[1].split('or')[1].replace(')', '').strip()}.",
                        html.Br(),
                        html.Span(html.Small(f"Can be observed with the {arrays_with_band(default_arrays, a_band)}."),
                            style={'color': '#888888'})
                    ], className="card-text"),
                ]), className="col-sm-3 my-2 shadow-1-strong")
            ]


# @app.callback([Output('initial-timeselection-div-guess', 'hidden'),
#               Output('initial-timeselection-div-epoch', 'hidden')],
#               [Input('initial-timeselection', 'value')])
# def type_initial_time_selection(time_selection_selected):
#     """Modifies the hidden message related to the two options about how to pick the observing time.
#     """
#     print('checked')
#     return [time_selection_selected, not time_selection_selected]


@app.callback([Output('initial-timeselection-div-guess', 'hidden'),
              Output('initial-timeselection-div-epoch', 'hidden')],
              [Input('initial-timeselection', 'value')])
def type_initial_time_selection(time_selection_selected):
    """Modifies the hidden message related to the two options about how to pick the observing time.
    """
    return [time_selection_selected, not time_selection_selected]



@app.callback([Output('timeselection-div-guess', 'hidden'),
              Output('timeselection-div-epoch', 'hidden')],
              [Input('timeselection', 'value')])
def type_time_selection(time_selection_selected):
    """Modifies the hidden message related to the two options about how to pick the observing time.
    """
    return [time_selection_selected, not time_selection_selected]



@app.callback(Output('band', 'value'),
              Input('initial-band', 'value'), prevent_initial_call=True)
def band_from_initial(initial_value):
    return tuple(fs.bands)[initial_value] if initial_value is not None else dash.no_update


@app.callback(Output('array', 'value'),
              [Input(f'network-{network.lower()}', 'value') for network in default_arrays if network != 'e-EVN'],
              prevent_initial_call=True)
def array_from_initial(*selected_networks):
    return [network for (network,selected) in \
            zip([n for n in default_arrays if n != 'e-EVN'], selected_networks) if selected]


@app.callback(Output('e-EVN', 'value'),
              Input('initial-e-EVN', 'value'), prevent_initial_call=True)
def e_EVN_from_initial(initial_value):
    return initial_value if initial_value is not None else dash.no_update


@app.callback(Output('timeselection', 'value'),
              Input('initial-timeselection', 'value'), prevent_initial_call=True)
def timeselection_from_initial(initial_value):
    return initial_value if initial_value is not None else dash.no_update


@app.callback(Output('starttime', 'date'),
              Input('initial-starttime', 'date'), prevent_initial_call=True)
def starttime_from_initial(initial_value):
    return initial_value if initial_value is not None else dash.no_update


@app.callback(Output('starthour', 'value'),
              Input('initial-starthour', 'value'), prevent_initial_call=True)
def starthour_from_initial(initial_value):
    return initial_value if initial_value is not None else dash.no_update


@app.callback(Output('duration', 'value'),
              Input('initial-duration', 'value'), prevent_initial_call=True)
def duration_from_initial(initial_value):
    return initial_value if initial_value is not None else dash.no_update


@app.callback(Output('source', 'value'),
              Input('initial-source', 'value'), prevent_initial_call=True)
def source_from_initial(initial_value):
    return initial_value if initial_value is not None else dash.no_update


@app.callback([Output('subbands', 'value'),
               Output('channels', 'value'),
               Output('pols', 'value'),
               Output('inttime', 'value')],
              Input('is_line', 'value'), prevent_initial_call=True)
def line_cont_setup(is_line_exp):
    if is_line_exp is None:
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update

    if is_line_exp:
        return 1, 4096, 4, 2
    else:
        return 8, 32, 4, 2


@app.callback([Output('button-picknetwork', 'disabled'),
               Output('button-picknetwork', 'children')],
        [Input(f"network-{network.lower()}", 'value') for network in default_arrays if network != 'e-EVN'])
def continue_from_networks(*networks):
    """Verifies that the user has selected at least one VLBI network during the wizard screen
    """
    for n in networks:
        if True in n:
            return False, 'Continue'

    return True, 'Select network(s) to continue'
    # len([True for n in networks if True in n]) > 0 % Slower


@app.callback([Output('button-pickband', 'disabled'),
               Output('button-pickband', 'children')],
              Input('initial-band', 'value'),
              [State(f"network-{network.lower()}", 'value') for network in default_arrays if network != 'e-EVN'])
def continue_from_band(selected_band, *networks):
    """Verifies that the selected band can be observed by the given network.
    """
    for n,nname in zip(networks, [network for network in default_arrays if network != 'e-EVN']):
        if True in n:
            if tuple(fs.bands.keys())[selected_band] in default_arrays[nname]['observing_bands']:
                return False, 'Continue'

    return True, 'The selected network cannot observe at this band'


@app.callback([Output('button-picktimes', 'disabled'),
         Output('button-picktimes', 'children')],
        [Input('initial-timeselection', 'value'),
         Input('initial-starttime', 'date'),
         Input('initial-starthour', 'value'),
         Input('initial-duration', 'value'),
         Input('initial-source', 'value')])
def continue_from_times(time_selection, time_date, time_hour, time_duration, source):
    """Verifies that the user has selected and introduced the required data before continue.
    """
    if (source is None) or (not verify_recognized_source(source)):
        return True, 'Specify epoch and target before continue'

    if time_selection:
        if (time_date is not None) and (time_hour is not None) and (time_duration is not None):
            try:
                dummy = float(time_duration)
                return (True, 'Specify epoch and target before continue') if (dummy <= 0) or (dummy > 4*24) \
                        else (False, 'Continue')
            except:
                return True, 'Specify epoch and target before continue'
        else:
            return True, 'Specify epoch and target before continue'

    return False, 'Continue'



@app.callback(Output('main-window2', 'children'),
              Output('is_line', 'value'),
              [Input('button-pickband', 'n_clicks'),
               Input('button-picknetwork', 'n_clicks'),
               Input('button-picktimes', 'n_clicks'),
               Input('button-mode-continuum', 'n_clicks'),
               Input('button-mode-line', 'n_clicks')])
def intro_choices(clicks_pickband, clicks_picknetwork, clicks_picktimes, clicks_continuum, clicks_line):
    if clicks_picknetwork is not None:
        return choice_page('band'), dash.no_update
    elif clicks_pickband is not None:
        return choice_page('time'), dash.no_update
    elif clicks_picktimes is not None:
        #TODO: if dates are provided, all date/time/dur must be provided
        return choice_page('mode'), dash.no_update
    elif clicks_continuum is not None:
        return choice_page('final'), False
    elif clicks_line is not None:
        return choice_page('final'), True
    else:
        return dash.no_update, dash.no_update



def initial_page():
    """Initial window with the two options to select: guided or manual setup of the observation.
    """
    return [
         html.Div(className='row justify-content-center', id='main-window2',
            children=html.Div(className='col-sm-6 justify-content-center',
                children=[html.Div(className='justify-content-center',
                    children=[#html.H3("Welcome!"),
                        html.P(["The EVN Observation Planner allows you to plan observations with the ",
                        html.A(href="https://www.evlbi.org", children="European VLBI Network"),
                        " (EVN) and other Very Long Baseline Interferometry (VLBI) networks. "
                        "The EVN Observation Planner helps you to determine when your source "
                        "can be observed by the different antennas, and provides the expected "
                        "outcome of these observations, like the expected sensitivity or resolution."]),
                        html.Br(),
                        html.Div(ge.initial_window_start(app))
                    ])
                ])
        )]


@app.callback(Output('full-window', 'children'),
              [Input('button-initial-wizard', 'n_clicks'),
               Input('button-initial-expert', 'n_clicks')])
def choice_for_setup(do_wizard, do_expert):
    if (do_expert is not None) or (do_wizard is not None):
        return [
            # order inverted to improve loading times
            html.Div(id='main-window2', hidden=do_expert is not None,
                     children=[dbc.Checklist(id='is_line', options=[{'label': 'line obs', 'value': False}],
                               value=[])] if do_expert is not None else choice_page('network')),
            html.Div(id='main-window', hidden=do_expert is None,
                     children=main_page(show_compute_button=do_expert is not None))
            ]
    else:
        return dash.no_update



def choice_page(choice_card):
    """Initial window with the introduction to the EVN Observation Planner and the band selection.
    """
    return [
        html.Div(className='row justify-content-center', id='main-window2',
            children=html.Div(className='col-sm-6 justify-content-center',
                children=[html.Div(className='justify-content-center',
                        children=[#html.H3("Welcome!"),
                                  html.P(["The EVN Observation Planner allows you to plan observations with the ",
                            html.A(href="https://www.evlbi.org", children="European VLBI Network"),
                            " (EVN) and other Very Long Baseline Interferometry (VLBI) networks. "
                            "The EVN Observation Planner helps you to determine when your source "
                            "can be observed by the different antennas, and provides the expected "
                            "outcome of these observations, like the expected sensitivity or resolution."]),
                            html.Br(),
                            *[
                                # html.Div(hidden=False if choice_card == 'choice' else True,
                                #          children=ge.initial_window_start(app)),
                                html.Div(hidden=False if choice_card == 'band' else True,
                                         children=ge.initial_window_pick_band()),
                                html.Div(hidden=False if choice_card == 'network' else True,
                                         children=ge.initial_window_pick_network(app, default_arrays)),
                                html.Div(hidden=False if choice_card == 'time' else True,
                                         children=ge.initial_window_pick_time()),
                                html.Div(hidden=False if choice_card == 'mode' else True,
                                         children=ge.initial_window_pick_mode(app)),
                                html.Div(hidden=False if choice_card == 'final' else True,
                                         children=ge.initial_window_final()),
                            ],
                        ], style={'text:align': 'justify !important'})
                ])
            )]




def main_page(results_visible=False, summary_output=None, fig_elev_output=None,
              fig_ant_output=None, fig_uv_output=None, fig_dirty_map_output=False, show_compute_button=True):
    return [# First row containing all buttons/options, list of telescopes, and button with text output
    dcc.ConfirmDialog(id='global-error', message=''),
    # Elements in second column (checkboxes with all stations)
    html.Div(className='container-fluid', children=[
        html.Div(className='row justify-content-center', children=[
             html.Div(className='col-sm-3', style={'max-width': '350px','float': 'left',
                                   'min-width': '17rem'}, children=[
                html.Div(className='form-group', children=[
                    html.Div('', style={'height': '70px'}),
                    dcc.Loading(id="loading2", children=[html.Div(id="loading-output2"), html.Br()],
                                type="dot"),
                    html.Button('Compute observation',
                        id='antenna-selection-button',
                        className='btn btn-primary btn-lg',
                        style={'width': '100%', 'margin-bottom': '1rem'})#if show_compute_button else html.Br(),
                ]),
                html.Br(),
                html.Div(className='form-group', children=[
                    html.H6(['Your observing Band',
                        *ge.tooltip(idname='popover-band',
                            message="This will update the "
                                    "antenna list showing the ones that can observe "
                                    "at that given frequency.")
                    ]),
                    dcc.Dropdown(id='band', persistence=True, value='18cm',
                         options=[{'label': fs.bands[b], 'value': b} for b \
                        # in fs.bands], value='18cm'),
                        in fs.bands], placeholder='Select observing band...'),
                ]),
                html.Br(),
                html.Div(className='form-group', children=[
                    html.H6(['Real-time correlation?',
                        *ge.tooltip(idname='popover-eevn',
                        message="Only available for the EVN: real-time correlation mode."
                                "The data are transferred and correlated in real-time, but "
                                "not all telescopes are capable for this and the bandwidth "
                                "may be limited. Observations during the e-EVN epochs.")
                    ]),
                    dbc.Checklist(id='e-EVN', className='checkbox', persistence=True,
                                  options=[{'label': ' e-EVN mode',
                                            'value': 'e-EVN'}], value=[]),
                ]),
                html.Br(),
                html.Div(className='form-group', children=[
                    html.H6(['Source (name or coordinates)',
                        *ge.tooltip(idname='popover-target',
                             message="Source name or coordinates. " \
                                     "You may see an error if the given name is not properly resolved. "
                                     "J2000 coordinates are assumed in both forms: 00:00:00 00:00:00 or " \
                                     "00h00m00s 00d00m00s.")
                    ]),
                    dcc.Input(id='source', value=None, type='text',
                              className='form-control', placeholder="hh:mm:ss dd:mm:ss",
                              persistence=True),
                    html.Small(id='error_source',
                               className='form-text text-muted'),
                ]),
                html.Br(),
                html.Div(className='form-group', children=[
                    html.H6('Epoch for observation'),
                    dbc.FormGroup([
                        dbc.RadioItems(options=[{"label": "I don't have a preferred epoch", "value": False},
                                                {"label": "I know the observing epoch", "value": True}],
                                       value=True, id="timeselection", inline=True, persistence=True),
                    ], inline=True),
                    html.Div(children=[
                        html.Div(id='timeselection-div-guess', className='row justify-content-center',
                            hidden=False, children=[
                                html.Small("Choose the first option to find out when your source "
                                           "may be visible (by >3 telescopes).", style={'color': '#999999'}),
                                html.Small("Note that this option may not provide the best (expected) "
                                           "results in case of combining different networks very far apart "
                                           "(e.g. LBA and EVN).", style={'color': '#999999'})
                        ]),
                        html.Div(id='timeselection-div-epoch', hidden=True, children=[
                            html.Label('Start of observation (UTC)'),
                            *ge.tooltip(idname='popover-startime', message="Select the date and "
                                        "time of the start of the observation (Universal, UTC, "
                                        "time). You will also see the day of the year (DOY) in "
                                        "brackets once the date is selected."),
                            html.Br(),
                            dcc.DatePickerSingle(id='starttime', date=None, min_date_allowed=dt(1900, 1, 1),
                                                 max_date_allowed=dt(2100, 1, 1),
                                                 display_format='DD-MM-YYYY (DDD)',
                                                 placeholder='Start date',
                                                 first_day_of_week=1,
                                                 initial_visible_month=dt.today(),
                                                 persistence=True,
                                                 className='form-picker'),
                            dcc.Dropdown(id='starthour', placeholder="Start time (UTC)", value=None,
                                         options=[{'label': f"{hm//60:02n}:{hm % 60:02n}", \
                                                   'value': f"{hm//60:02n}:{hm % 60:02n}"} \
                                                  for hm in range(0, 24*60, 15)],
                                         persistence=True, className='form-hour'),
                            html.Small(id='error_starttime', style={'color': 'red'},
                                       className='form-text text-muted'),
                            html.Label('Duration of the observation (in hours)'),
                            html.Div(className='form-group', children=[
                                dcc.Input(id='duration', value=None, type='number', className='form-control',
                                           placeholder="Duration in hours", persistence=True, inputMode='numeric'),
                                html.Small(id='error_duration', className='form-text text-danger')
                            ])
                        ])
                    ]),
                ]),
                html.Br(),
                html.Div(className='form-group', children=[
                    html.H6(['% of on-target time',
                        *ge.tooltip(idname='popover-ontarget',
                             message="Assumes that you will only spend this amount of the total " \
                                     "observing time on the given target source. It affects the " \
                                     "expected sensitivity."),
                    ]),
                    dcc.Slider(id='onsourcetime', min=20, max=100, step=5, value=70,
                               marks= {i: str(i) for i in range(20, 101, 10)},
                               persistence=True),
                    html.Label(id='onsourcetime-label', style={'color': '#999999'},
                               children='70% of the observation.'),
                ]),
                html.Br(),
                html.Div(className='form-group', children=[
                    html.H6(['Datarate per station',
                        *ge.tooltip(idname='popover-datarate',
                             message=["Expected datarate for each station, assuming all " \
                                      "of them run at the same rate.",
                                 html.Ul([
                                    html.Li("The EVN can run typically at up to 2 Gbps (1 Gbps at L band), " \
                                            "although a few antennas may observe at lower datarates."),
                                    html.Li("The VLBA can now observe up to 4 Gbps."),
                                    html.Li("The LBA typically runs at 512 Mbps but can reach up to 1 Gbps."),
                                    html.Li("Check the documentation from other networks to be " \
                                            "sure about their capabilities.")])])
                    ]),
                    dcc.Dropdown(id='datarate',
                                 placeholder="Select the data rate...",
                                 options=[{'label': fs.data_rates[dr], 'value': dr} \
                                 for dr in fs.data_rates], value=2048, persistence=True),
                    html.Label(id='bandwidth-label', style={'color': '#999999'}, children='')
                ]),
                html.Div(className='form-group', children=[
                    html.H6(['Number of subbands',
                        *ge.tooltip(idname='popover-subbands',
                             message="Number of subbands to split the total observed bandwidth "
                                     " during correlation (IFs in AIPS).")
                    ]),
                    dcc.Dropdown(id='subbands', placeholder="Select no. subbands...",
                                 options=[{'label': fs.subbands[sb], 'value': sb} \
                                 for sb in fs.subbands], value=8, persistence=True),
                ]),
                html.Br(),
                html.Div(className='form-group', children=[
                    html.H6(['Number of spectral channels',
                        *ge.tooltip(idname='popover-channels',
                             message="How many channels per subband will be produced "
                                     "during correlation.")
                    ]),
                    dcc.Dropdown(id='channels', placeholder="Select no. channels...",
                                 options=[{'label': fs.channels[ch],
                                           'value': ch} \
                                 for ch in fs.channels], value=32, persistence=True),
                ]),
                html.Br(),
                html.Div(className='form-group', children=[
                    html.H6(['Number of polarizations',
                        *ge.tooltip(idname='popover-pols',
                        message="Number of polarizations to correlate. Note that VLBI uses circular " \
                                "polarizations. Full polarization implies the four stokes: RR, LL, RL, LR; " \
                                "while dual polarization implies RR and LL only.")
                    ]),
                    dcc.Dropdown(id='pols', placeholder="Select polarizations...",
                                 options=[{'label': fs.polarizations[p], 'value': p} \
                                 for p in fs.polarizations], value=4, persistence=True),
                ]),
                html.Br(),
                html.Div(className='form-group', children=[
                    html.H6(['Integration time',
                        *ge.tooltip(idname='popover-inttime',
                        message="Integration time to compute each visibility. Note that for continuum " \
                                "observations values of 1-2 seconds are typical.")
                    ]),
                    dcc.Dropdown(id='inttime', placeholder="Select integration time...",
                                 options=[{'label': fs.inttimes[it], 'value': it} \
                                 for it in fs.inttimes], value=2, persistence=True),
                ])
            ]),
            dcc.Tabs(parent_className='custom-tabs col', className='custom-tabs-container', id='tabs',
                value='tab-setup', children=[
                dcc.Tab(label='Antenna Selection', className='custom-tab', value='tab-setup',
                        selected_className='custom-tab--selected', children=[
                    # Elements in first column
                    html.Div(className='row justify-content-center', children=[
                    html.Div(className='col-9', children=[
                        html.Div(id='first-advise', className='col-sm-9', children=[
                            html.H4("Customize your observation"),
                            html.P(["Select which VLBI network(s) you want to use in your "
                                   "observation, or select an ad-hoc array of antennas. ", html.Br(),
                                   "Set the basic information from your observations: "
                                   "observing band, target source, epoch, and observing setup. ", html.Br(),
                                   "Finally, press the blue ", html.B("'compute observation'"),
                                   " button. ", html.Br(),
                                   "You will get a detailed "
                                   "summary of the planned observation and expected outcome in the different "
                                   "tabs."]),
                            html.P(html.Em(["Note that only antennas that can observe at the selected band "
                                    "will be clickable."], className='form-text text-warning'))
                        ], style={'margin-top': '2rem', 'margin-bottom': '2rem'}),
                        html.Div(className='col-9 form-group row align-items-end', children=[
                            html.Div(className='col-md-12', children=[
                                html.Label('Select default VLBI Network(s)'),
                                        # style={'color': '#a01d26'}),
                                *ge.tooltip(idname='popover-network', message="Automatically selects "
                                            "the default participating antennas for the selected VLBI network(s)."),
                                dcc.Dropdown(id='array', options=[{'label': f"{n}: {default_arrays[n]['name']}" \
                                        if n != default_arrays[n]['name'] else n,
                                        'value': n} for n in default_arrays if n != 'e-EVN'], value=[],
                                        persistence=True, multi=True),
                            ]),
                        ]),
                        html.Div(className='col-9 text-center justify-content-center', children=[
                            dcc.Loading(id="loading", children=[html.Br(), html.Div(id="loading-output")],
                                        type="dot")
                        ]),
                        html.Div([dbc.Tooltip(ge.antenna_card(app, s), placement='right',
                            hide_arrow=True, target=f"_input_{s.codename}",
                            innerClassName='tooltip-card-inner') for s in all_antennas
                        ]),
                        html.Div(id='antennas-div', className='container', children=[
                            html.Div(className='antcheck', children=[html.Br(), html.Br(),
                                html.Label(html.H4(f"{sorted_networks[an_array]}"),
                                           style={'width': '100%'}),
                                html.Br(),
                                html.Div(className='antcheck', children=[
                                    dbc.FormGroup([
                                        Checkbox(id=f"check_{s.codename}", persistence=True,
                                                 className='custom-control-input',
                                                 disabled=not s.has_band(selected_band)),
                                        dbc.Label(s.name, html_for=f"check_{s.codename}",
                                                  id=f"_input_{s.codename}",
                                                 className='custom-control-label')
                                    ], check=True, inline=True,
                                    className="custom-checkbox custom-control custom-control-inline")
                                    for s in all_antennas if s.network == an_array
                                ])
                            ]) for an_array in sorted_networks
                        ]),
                        html.Div(style={'height': '15rem'})
                    ]),
                    # html.Div(className='col-sm-2', style={'float': 'left'}, children=[
                    # ])
                    ])
                ]),

                dcc.Tab(label='Summary', className='custom-tab', value='tab-summary', id='tab-summary',
                        selected_className='custom-tab--selected', disabled=not results_visible, children=[
                    html.Div(className='row justify-content-center', children=[
                        html.Div(className='col-10 justify-content-center',
                                 id='sensitivity-output',
                                 children=summary_output)
                                #  children=[html.Div(className='col-md-6', children=[
                                #     html.Br(), html.Br(), html.H2("Set the observation first"),
                                #     html.P("Here you will see a summary of your observation, "
                                #            "with information about all participating stations, longest and "
                                #            "shortest baseline, expected size of the data once is correlated, "
                                #            "reached resolution and sensitivity, and the limitations in your "
                                #            "field of view due to time and frequency smearing.")])
                                # ])
                    ])
                ]),
                dcc.Tab(label='Elevations', className='custom-tab', value='tab-elevation', id='tab-elevation',
                        selected_className='custom-tab--selected',
                        disabled=not results_visible, children=[
                    html.Div(className='row justify-content-center', children=[
                    html.Div(className='col-md-8 justify-content-center', children=[
                        # Elevation VS time
                        html.Br(),
                        html.Div([
                            html.Br(),
                            html.H4("When is your source visible?"),
                            html.Br(),
                            dbc.Alert([html.H4("Interactive plots", className='alert-heading'),
                                       html.P(["A single click on one station in the legend will "
                                               "hide/show it.", html.Br(), "Double-click will hide/show "
                                               "all other antennas. You can also save the plot "
                                               "as png."]),
                                      ], color='info', dismissable=True),
                            html.Br(),
                            html.P("The following plot shows the source elevation for the "
                            "different antennas during the proposed observation. The horizontal "
                            "solid and dashed lines represent the elevation of 20 and 10 degrees, "
                            "respectively.")
                        ]),
                        html.Div([
                            dcc.Graph(id='fig-elev-time', children=fig_elev_output)
                        ],className='tex2jax_ignore'),
                        html.Br(),
                        html.Div([
                            html.P("""The following plot shows when the source may be observed
                            for the different antennas, assuming a minimum elevation of 10 degrees
                            for most antennas (except e.g. Arecibo). Note that some antennas may
                            have additional constraints for particular azimuth or elevation
                            angles that are not considered here.
                            """)
                        ]),
                        html.Div([
                            dcc.Graph(id='fig-ant-time', children=fig_ant_output)
                        ],className='tex2jax_ignore')
                    ])])
                ]),
                dcc.Tab(label='UV Coverage', className='custom-tab', value='tab-uv', id='tab-uv',
                        selected_className='custom-tab--selected', disabled=fig_uv_output is None, children=[
                    #  Images
                    html.Div(className='row justify-content-center', children=[
                    html.Div(className='col-md-8 justify-content-center', children=[
                        # dcc.Markdown(children="""To be implemented.
                        #     The uv coverage and expected dirty images will go here.""")
                        html.Br(),
                        html.Div([
                            html.Br(),
                            html.H4("Resulting (u,v) coverage"),
                            html.Br()
                        ]),
                        html.Div(children=[dcc.Graph(id='fig-uvplane', children=fig_uv_output)],
                                 className='tex2jax_ignore'),
                        html.Br(),
                        html.Div([
                            html.Br(),
                            html.H4("Resulting dirty beam"),
                            html.Br()
                        ]),
                        html.Div(children=[dcc.Graph(id='fig-dirtymap', figure=fig_dirty_map_output)],
                                 className='tex2jax_ignore')
                    ])])
                ]),
                dcc.Tab(label='Documentation', className='custom-tab', value='tab-doc',
                        selected_className='custom-tab--selected', children=[
                    #  Documentation
                    html.Div(className='row justify-content-center', children=[
                        html.Div([html.Br(), html.Br()]),
                        html.Div(className='col-md-8', children=get_doc_text())
                    ])
                ])
            ])
        ]),
    html.Div(className='container-fluid', children=[html.Br(), html.Br()])
    ])
    ]



@app.callback(Output('onsourcetime-label', 'children'),
              [Input('onsourcetime', 'value'),
               Input('duration', 'value'),
               Input('timeselection', 'value')])
def update_onsourcetime_label(onsourcetime, total_duration, defined_epoch):
    """Keeps the on-source time label updated with the value selected by the user.
    """
    if (total_duration is not None) and defined_epoch:
        return f"{onsourcetime}% of the total time  ({ge.optimal_units(total_duration*u.h*onsourcetime/100, [u.h, u.min, u.s]):.03n})."

    return f"{onsourcetime}% of the total observation."



@app.callback(Output('bandwidth-label', 'children'),
        [Input('datarate', 'value'),
         Input('pols', 'value')])
def update_bandwidth_label(datarate, npols):
    """Updates the total bandwidth label as a function of the selected datarate and number of
    polarizations. Returns a string with the value and units.
    """
    if (None not in (datarate, npols)) and (datarate != -1):
        # Either 1 or 2 pols per station:
        temp = npols % 3 + npols // 3
        return [f"The total bandwidth is {ge.optimal_units(datarate*u.MHz/(temp*2*2), [u.GHz, u.MHz, u.kHz] )}.",
                html.Br(), html.Br()]

    return ''





@app.callback([Output(f"check_{s.codename}", 'checked') for s in all_antennas] + \
              [Output(f"check_{s.codename}", 'disabled') for s in all_antennas] + \
              [Output('datarate', 'value')],
              [Input('band', 'value'),
               Input('array', 'value'),
               Input('e-EVN', 'value'),
               Input('is_line', 'value')])
def select_antennas(selected_band, selected_networks, is_eEVN, is_line):
    """Given a selected band and selected default networks, it selects the associated
    antennas from the antenna list.
    """
    # Getting the data rate
    if 'EVN' in selected_networks or is_eEVN:
        datarate = default_arrays['EVN']['max_datarate'] if selected_band not in ('18cm', '21cm') else 1024
    else:
        datarate = -1
        for an_array in selected_networks:
            datarate = max(datarate, default_arrays[an_array]['max_datarate'])

    # Getting the selected antennas
    selected_antennas = []
    if is_eEVN:
        selected_antennas = [ant for ant in default_arrays['e-EVN']['default_antennas'] \
                         if all_antennas[ant].has_band(selected_band)]

    for an_array in selected_networks:
            selected_antennas += [ant for ant in default_arrays[an_array]['default_antennas'] \
                if all_antennas[ant].has_band(selected_band) and (an_array != 'EVN' or not is_eEVN)]

    return [True if s.codename in selected_antennas else False for s in all_antennas] + \
           [False if (s.has_band(selected_band) and (s.real_time or not is_eEVN)) else True \
            for s in all_antennas] + [datarate if not is_line else 256]



@app.callback([Output('initial-error_starttime', 'children'),
               Output('initial-error_duration', 'children')],
              [Input('initial-starttime', 'date'), Input('starthour', 'value'),
               Input('initial-duration', 'value')])
def check_initial_obstime(starttime, starthour, duration):
    """Verify the introduced times/dates for correct values.
    Once the user has introduced all values for the start and end of the observation,
    it guarantees that they have the correct shape:
        - the duration of the observation is > 0 hours.
        - The total observing length is less than five days (value chosen for computational reasons).
    """
    if duration is None:
        return "", ""

    if (not isinstance(duration, float)) and (not isinstance(duration, int)):
        return "", "Must be a number"

    if duration <= 0.0:
        return "", "The duration must be a positive number"
    elif duration > 4*24:
        return "", "Please, put an observation shorter than 4 days"

    return "", ""


@app.callback([Output('error_starttime', 'children'),
               Output('error_duration', 'children')],
              [Input('starttime', 'date'), Input('starthour', 'value'),
               Input('duration', 'value')])
def check_obstime(starttime, starthour, duration):
    """Verify the introduced times/dates for correct values.
    Once the user has introduced all values for the start and end of the observation,
    it guarantees that they have the correct shape:
        - the duration of the observation is > 0 hours.
        - The total observing length is less than five days (value chosen for computational reasons).
    """
    if duration is None:
        return "", ""

    if (not isinstance(duration, float)) and (not isinstance(duration, int)):
        return "", "Must be a number"

    if duration <= 0.0:
        return "", "The duration must be a positive number"
    elif duration > 4*24:
        return "", "Please, put an observation shorter than 4 days"

    return "", ""


@app.callback([Output('initial-error_source', 'children'),
        Output('initial-error_source', 'className')],
        [Input('initial-source', 'value')])
def get_initial_source(source_coord):
    """Verifies that the introduced source coordinates have a right format.
    If they are correct, it does nothing. If they are incorrect, it shows an error label.
    """

    if source_coord != 'hh:mm:ss dd:mm:ss' and source_coord != None and source_coord != '':
        if len(source_coord) > 30:
            # Otherwise the source name check gets too slow
            return "Name too long.", 'form-text text-danger'
        try:
            dummy_target = observation.Source(convert_colon_coord(source_coord), 'Source')
            return '', dash.no_update
        except ValueError as e:
            try:
                dummy_target = coord.get_icrs_coordinates(source_coord)
                return dummy_target.to_string('hmsdms'), 'form-text text-muted'
            except coord.name_resolve.NameResolveError as e:
                return "Unrecognized name. Use 'hh:mm:ss dd:mm:ss' or 'XXhXXmXXs XXdXXmXXs'", \
                       'form-text text-danger'
    else:
        return '', dash.no_update


def verify_recognized_source(a_source):
    """Equivalent to the previous function, but returns a bool for when a source name/coordinates
    have been introduced correctly or not.
    """
    if a_source is None:
        return False

    if len(a_source) > 30:
        return False

    try:
        dummy_target = observation.Source(convert_colon_coord(a_source), 'Source')
        return True
    except ValueError as e:
        try:
            dummy_target = coord.get_icrs_coordinates(a_source)
            return True
        except coord.name_resolve.NameResolveError as e:
            return False

    return False



@app.callback([Output('error_source', 'children'),
        Output('error_source', 'className')],
        [Input('source', 'value')])
def get_source(source_coord):
    """Verifies that the introduced source coordinates have a right format.
    If they are correct, it does nothing. If they are incorrect, it shows an error label.
    """

    if source_coord != 'hh:mm:ss dd:mm:ss' and source_coord != None and source_coord != '':
        if len(source_coord) > 30:
            # Otherwise the source name check gets too slow
            return "Name too long.", 'form-text text-danger'
        try:
            dummy_target = observation.Source(convert_colon_coord(source_coord), 'Source')
            return '', dash.no_update
        except ValueError as e:
            try:
                dummy_target = coord.get_icrs_coordinates(source_coord)
                return dummy_target.to_string('hmsdms'), 'form-text text-muted'
            except coord.name_resolve.NameResolveError as e:
                return "Unrecognized name. Use 'hh:mm:ss dd:mm:ss' or 'XXhXXmXXs XXdXXmXXs'", \
                       'form-text text-danger'
    else:
        return '', dash.no_update



@app.callback([Output('loading-output', 'children'),
               Output('loading-output2', 'children'),
               Output('main-window', 'hidden'),
               Output('main-window2', 'hidden'),
               Output('sensitivity-output', 'children'),
               Output('fig-elev-time', 'figure'),
               Output('fig-ant-time', 'figure'),
               Output('fig-uvplane', 'figure'),
               Output('fig-dirtymap', 'figure'),
               Output('global-error', 'message'),
               Output('tabs', 'value'),
               Output('tab-summary', 'disabled'),
               Output('tab-elevation', 'disabled'),
               Output('tab-uv', 'disabled')],
              [Input('antenna-selection-button', 'n_clicks')],
              [State('band', 'value'),
               State('starttime', 'date'),
               State('starthour', 'value'),
               State('duration', 'value'),
               State('source', 'value'),
               State('onsourcetime', 'value'),
               State('datarate', 'value'),
               State('subbands', 'value'),
               State('channels', 'value'),
               State('pols', 'value'),
               State('inttime', 'value'),
               State('timeselection', 'value'),
               State('tabs', 'value')] + \
               [State(f"check_{s.codename}", 'checked') for s in all_antennas])
def compute_observation(n_clicks, band, starttime, starthour, duration, source, onsourcetime,
                        datarate, subbands, channels, pols, inttime, epoch_selected, selected_tab, *ants):
    """Computes all products to be shown concerning the set observation.
    """
    # To decide where to put the output message
    out_center = selected_tab == 'tab-setup' or selected_tab == 'tab-doc'
    if n_clicks is None:
        return '', '', dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, \
               dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, \
               dash.no_update, dash.no_update

    if epoch_selected:
        # All options must be completed
        if None in (band, starttime, starthour, duration, source, datarate, subbands, channels, pols, inttime) \
                or source == "":
            missing = [label for label,atr in zip(('observing band', 'target source', 'start observing date',
                        'start observing time', 'duration of the observation', 'data rate', 'number of subbands',
                        'number of channels', 'number of polarizations', 'integration time'), (band, source,
                        starttime, starthour, duration, datarate, subbands, channels, pols, inttime)) \
                        if (atr is None) or (atr == "")]
            temp = [alert_message(["Complete all fields and options before computing the observation.\n" + \
                                  f"Currently it is missing: {', '.join(missing)}."]), '']
            return *[temp if out_center else temp[::-1]][0], dash.no_update, dash.no_update, dash.no_update, \
                    dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, \
                    dash.no_update, dash.no_update, dash.no_update, dash.no_update
    else:
        # All options but the ones related to the observing epoch must be completed
        if None in (band, source, datarate, subbands, channels, pols, inttime) or source == "":
            missing = [label for label,atr in zip(('observing band', 'target source', 'data rate',
                        'number of subbands', 'number of channels', 'number of polarizations',
                        'integration time'), (band, source, starttime, starthour, duration, datarate, subbands,
                        channels, pols, inttime)) if (atr is None) or (atr == "")]
            temp = [alert_message(["Complete all fields and options before computing the observation.\n" + \
                                  f"Currently it is missing: {', '.join(missing)}."]), '']
            return *[temp if out_center else temp[::-1]][0], dash.no_update, dash.no_update, dash.no_update, \
                    dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, \
                    dash.no_update, dash.no_update, dash.no_update, dash.no_update

    if ants.count(True) == 0:
        temp = [alert_message(["You need to select the antennas you wish to observe your source. " \
                "Either manually or by selected a default VLBI network at your top left."]), '']
        return *[temp if out_center else temp[::-1]][0], dash.no_update, dash.no_update, dash.no_update, \
                dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, \
                dash.no_update, dash.no_update, dash.no_update, dash.no_update
    # A single antenna computation is not supported
    if ants.count(True) == 1:
        temp = [alert_message(["Single-antenna computations are not suported. " \
                              "Please choose at least two antennas"]), dash.no_update]
        return *[temp if out_center else temp[::-1]][0], dash.no_update, dash.no_update, dash.no_update, \
                dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, \
                dash.no_update, dash.no_update, dash.no_update, dash.no_update

    try:
        target_source = observation.Source(convert_colon_coord(source), 'Source')
    except ValueError as e:
        try:
            target_source = observation.Source(coord.get_icrs_coordinates(source), source)
        except coord.name_resolve.NameResolveError as e:
            temp = [alert_message(["Wrong source name or coordinates.", html.Br(),
                    "Either the source name hasn't been found or the coordinates format is incorrect."]), \
                    "First, set correctly an observation in the previous tab.", '']
            return *[temp if out_center else temp[::-1]][0], dash.no_update, dash.no_update, dash.no_update, \
                    dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, \
                    dash.no_update, dash.no_update, dash.no_update, dash.no_update
    if not epoch_selected:
        try:
            utc_times, _ = observation.Observation.guest_times_for_source(target_source,
                            stations.Stations('Observation', itertools.compress(all_antennas, ants)))
        except observation.SourceNotVisible:
            temp = [alert_message([
                        html.P(["Your source cannot be observed within the arranged observation.",
                        html.Br(),
                        "There are no antennas that can simultaneously observe your source "
                        "during the given observing time."]),
                        html.P("Modify the observing time or change the selected antennas"
                               " to observe this source.")], title="Warning!"), '']
            return *[temp if out_center else temp[::-1]][0], dash.no_update, dash.no_update, dash.no_update, \
                    dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, \
                    dash.no_update, dash.no_update, dash.no_update, dash.no_update

        obs_times = utc_times[0] + np.linspace(0, (utc_times[1]-utc_times[0]).to(u.min).value, 50)*u.min
    else:
        try:
            time0 = Time(dt.strptime(f"{starttime} {starthour}", '%Y-%m-%d %H:%M'),
                         format='datetime', scale='utc')
        except ValueError as e:
            temp = [alert_message("Incorrect format for starttime."), \
                   "First, set correctly an observation in the previous tab.", '']
            return *[temp if out_center else temp[::-1]][0], dash.no_update, dash.no_update, dash.no_update, \
                    dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, \
                    dash.no_update, dash.no_update, dash.no_update, dash.no_update

        if duration <= 0.0:
            temp = [alert_message("The duration of the observation must be a positive number of hours"), \
                   "First, set correctly an observation in the previous tab.", '']
            return *[temp if out_center else temp[::-1]][0], dash.no_update, dash.no_update, dash.no_update, \
                    dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, \
                    dash.no_update, dash.no_update, dash.no_update, dash.no_update

        if duration > 4*24.0:
            temp = [alert_message("Please, set an observation that lasts for less than 4 days."), \
                   "First, set correctly an observation in the previous tab.", '']
            return *[temp if out_center else temp[::-1]][0], dash.no_update, dash.no_update, dash.no_update, \
                    dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, \
                    dash.no_update, dash.no_update, dash.no_update, dash.no_update

        obs_times = time0 + np.linspace(0, duration*60, 50)*u.min

    try:
        obs = observation.Observation(target=target_source, times=obs_times, band=band,
                      datarate=datarate, subbands=subbands, channels=channels,
                      polarizations=pols, inttime=inttime, ontarget=onsourcetime/100.0,
                      stations=stations.Stations('Observation',
                                                 itertools.compress(all_antennas, ants)))
        sensitivity_results = update_sensitivity(obs)
    except observation.SourceNotVisible:
        temp = [alert_message([
                    html.P(["Your source cannot be observed within the arranged observation.",
                    html.Br(),
                    "There are no antennas that can simultaneously observe your source "
                    "during the given observing time."]),
                    html.P("Modify the observing time or change the selected antennas"
                           " to observe this source.")], title="Warning!"), '']
        return *[temp if out_center else temp[::-1]][0], dash.no_update, dash.no_update, dash.no_update, \
                dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, \
                dash.no_update, dash.no_update, dash.no_update, dash.no_update

    # TODO: parallelize all these fig functions
    if out_center:
        return '', '', False, True, sensitivity_results, get_fig_ant_elev(obs), \
            get_fig_ant_up(obs), get_fig_uvplane(obs), get_fig_dirty_map(obs), dash.no_update, \
           'tab-summary', False, False, False
    else:
        return '', dbc.Alert("Results have been updated.", color='info', dismissable=True), False, True, \
           sensitivity_results, get_fig_ant_elev(obs), get_fig_ant_up(obs), get_fig_uvplane(obs), \
           get_fig_dirty_map(obs), dash.no_update, \
           dash.no_update, False, False, False



def get_fig_ant_elev(obs):
    data_fig = []
    data_dict = obs.elevations()
    # Some reference lines at low elevations
    for ant in data_dict:
        data_fig.append({'x': obs.times.datetime, 'y': data_dict[ant].value,
                        'mode': 'lines', 'hovertemplate': "Elev: %{y:.2n}º (%{x})",
                        'name': obs.stations[ant].name})

    data_fig.append({'x': obs.times.datetime, 'y': np.zeros_like(obs.times)+10,
                     'mode': 'lines', 'hoverinfo': 'skip', 'name': 'Elev. limit 10º',
                     'line': {'dash': 'dash', 'opacity': 0.5, 'color': 'gray'}})
    data_fig.append({'x': np.unwrap(obs.gstimes.value*2*np.pi/24)*24/(2*np.pi), 'y': np.zeros_like(obs.times)+20,
                     'xaxis': 'x2', 'mode': 'lines', 'hoverinfo': 'skip',
                     'name': 'Elev. limit 20º', 'line': {'dash': 'dot', 'opacity': 0.5,
                     'color': 'gray'}})
    return {'data': data_fig,
            'layout': {'title': 'Source elevation during the observation',
                       'hovermode': 'closest',
                       'xaxis': {'title': 'Time (UTC)', 'showgrid': False,
                                 'ticks': 'inside', 'showline': True, 'mirror': False,
                                 'hovermode': 'closest', 'color': 'black'},
                       'xaxis2': {'title': {'text': 'Time (GST)', 'standoff': 0},
                                  'showgrid': False, 'overlaying': 'x', #'dtick': 1.0,
                                  'tickvals': np.arange(np.ceil(obs.gstimes.value[0]),
                                        np.floor(np.unwrap(obs.gstimes.value*2*np.pi/24)[-1]*24/(2*np.pi))+1),
                                  'ticktext': np.arange(np.ceil(obs.gstimes.value[0]),
                                        np.floor(np.unwrap(obs.gstimes.value*2*np.pi/24)[-1]*24/(2*np.pi))+1)%24,
                                  'ticks': 'inside', 'showline': True, 'mirror': False,
                                  'hovermode': 'closest', 'color': 'black', 'side': 'top'},
                       'yaxis': {'title': 'Elevation (degrees)', 'range': [0., 92.],
                                 'ticks': 'inside', 'showline': True, 'mirror': "all",
                                 'showgrid': False, 'hovermode': 'closest'},
                       'zeroline': True, 'zerolinecolor': 'k'}}



def get_fig_ant_up(obs):
    data_fig = []
    data_dict = obs.is_visible()
    gstimes = np.unwrap(obs.gstimes.value*2*np.pi/24)*24/(2*np.pi)
    gstimes = np.array([dt(obs.times.datetime[0].year, obs.times.datetime[0].month, obs.times.datetime[0].day) \
                        + datetime.timedelta(seconds=gst*3600) for gst in gstimes])
    for i,ant in enumerate(data_dict):
        # xs = [obs.times.datetime[0].date() + datetime.timedelta(seconds=i*3600) for i in np.unwrap(obs.gstimes.value*2*np.pi/24)[data_dict[ant]]*24/(2*np.pi)]
        xs = gstimes[data_dict[ant]]
        data_fig.append({'x': xs,
                         'y': np.zeros_like(data_dict[ant][0])-i, 'type': 'scatter',
                         'hovertemplate': "GST %{x}",
                         'mode': 'markers', 'marker_symbol': "41",
                         'hoverinfo': "skip",
                         'name': obs.stations[ant].name})

    data_fig.append({'x': obs.times.datetime, 'y': np.zeros_like(obs.times)-0.5,
                     'xaxis': 'x2',
                     'mode': 'lines', 'hoverinfo': 'skip', 'showlegend': False,
                     'line': {'dash': 'dot', 'opacity': 0.0, 'color': 'white'}})
    return {'data': data_fig,
            'layout': {'title': {'text': 'Source visible during the observation',
                                 'y': 1, 'yanchor': 'top'},
                       'hovermode': 'closest',
                       'xaxis': {'title': 'Time (GST)', 'showgrid': False,
                                 'range': [gstimes[0], gstimes[-1]],
                                 # 'tickvals': np.arange(np.ceil(obs.gstimes.value[0]),
                                 #        np.floor(np.unwrap(obs.gstimes.value*2*np.pi/24)[-1]*24/(2*np.pi))+1),
                                 # 'ticktext': np.arange(np.ceil(obs.gstimes.value[0]),
                                 #        np.floor(np.unwrap(obs.gstimes.value*2*np.pi/24)[-1]*24/(2*np.pi))+1)%24,
                                 'tickformat': '%H:%M',
                                 'ticks': 'inside', 'showline': True, 'mirror': False,
                                 'hovermode': 'closest', 'color': 'black'},
                       'xaxis2': {'title': {'text': 'Time (UTC)', 'standoff': 0},
                                  'showgrid': False, 'overlaying': 'x', #'dtick': 1.0,
                                  'ticks': 'inside', 'showline': True, 'mirror': False,
                                  'hovermode': 'closest', 'color': 'black', 'side': 'top'},
                       'yaxis': {'ticks': '', 'showline': True, 'mirror': True,
                                 'showticklabels': False, 'zeroline': False,
                                 'showgrid': False, 'hovermode': 'closest',
                                 'startline': False}}}

    data_fig = []
    data_dict = obs.is_visible()
    for i,ant in enumerate(data_dict):
        data_fig.append({'x': obs.times.datetime[data_dict[ant]],
                         'y': np.zeros_like(data_dict[ant][0])-i, 'type': 'scatter',
                         'hovertemplate': "%{x}",
                         'mode': 'markers', 'marker_symbol': "41",
                         'hoverinfo': "skip",
                         'name': obs.stations[ant].name})

    data_fig.append({'x': np.unwrap(obs.gstimes.value), 'y': np.zeros_like(obs.times)-0.5,
                     'xaxis': 'x2',
                     'mode': 'lines', 'hoverinfo': 'skip', 'showlegend': False,
                     'line': {'dash': 'dot', 'opacity': 0.0, 'color': 'white'}})
    return {'data': data_fig,
            'layout': {'title': 'Source visible during the observation',
                       'xaxis': {'title': 'Time (UTC)', 'showgrid': False,
                                 'ticks': 'inside', 'showline': True, 'mirror': False,
                                 'hovermode': 'closest', 'color': 'black'},
                       'xaxis2': {'title': {'text': 'Time (GST)', 'standoff': 0},
                                  'showgrid': False, 'overlaying': 'x', #'dtick': 1.0,
                                  'tickvals': np.arange(np.ceil(obs.gstimes.value[0]),
                                            np.floor(np.unwrap(obs.gstimes.value)[-1])+1),
                                  'ticktext': np.arange(np.ceil(obs.gstimes.value[0]),
                                            np.floor(np.unwrap(obs.gstimes.value)[-1])+1) % 24,
                                  'ticks': 'inside', 'showline': True, 'mirror': False,
                                  'hovermode': 'closest', 'color': 'black', 'side': 'top'},
                       'yaxis': {'ticks': '', 'showline': True, 'mirror': True,
                                 'showticklabels': False, 'zeroline': False,
                                 'showgrid': False, 'hovermode': 'closest',
                                 'startline': False}}}


def get_fig_uvplane(obs):
    data_fig = []
    bl_uv = obs.get_uv_baseline()
    for bl_name in bl_uv:
        # accounting for complex conjugate
        uv = np.empty((2*len(bl_uv[bl_name]), 2))
        uv[:len(bl_uv[bl_name]), :] = bl_uv[bl_name]
        uv[len(bl_uv[bl_name]):, :] = -bl_uv[bl_name]
        data_fig.append({'x': uv[:,0],
                         'y': uv[:,1],
                         # 'type': 'scatter', 'mode': 'lines',
                         'type': 'scatter', 'mode': 'markers',
                         'marker': {'symbol': '.', 'size': 2},
                         'name': bl_name, 'hovertext': bl_name, 'hoverinfo': 'name', 'hovertemplate': ''})
    return {'data': data_fig,
            'layout': {'title': '', 'showlegend': False,
                       'hovermode': 'closest',
                       'width': 700, 'height': 700,
                       'xaxis': {'title': 'u (lambda)', 'showgrid': False, 'zeroline': False,
                                 'ticks': 'inside', 'showline': True, 'mirror': "all",
                                 'color': 'black'},
                       'yaxis': {'title': 'v (lambda)', 'showgrid': False, 'scaleanchor': 'x',
                                 'ticks': 'inside', 'showline': True, 'mirror': "all",
                                 'color': 'black', 'zeroline': False}}}



def get_fig_dirty_map(obs):
    dirty_map_nat, laxis = obs.get_dirtymap(pixsize=1024, robust='natural')
    dirty_map_uni, laxis = obs.get_dirtymap(pixsize=1024, robust='uniform')
    fig1 = px.imshow(img=dirty_map_nat, x=laxis, y=laxis[::-1], labels={'x': 'RA (mas)', 'y': 'Dec (mas)'}, \
            aspect='equal')
    fig2 = px.imshow(img=dirty_map_uni, x=laxis, y=laxis[::-1], labels={'x': 'RA (mas)', 'y': 'Dec (mas)'}, \
            aspect='equal')
    # for i,f in enumerate([fig1, fig2]):
    #     f.layout.xaxis.autorange = "reversed"
    #     f.layout.yaxis.autorange = True
    #     f['layout']['coloraxis']['showscale'] = False

    # fig1.layout.xaxis.autorange = "reversed"
    # fig2.layout.xaxis.autorange = "reversed"
    fig = make_subplots(rows=1, cols=2, subplot_titles=('Natural weighting', 'Uniform weighting'),
                        shared_xaxes=True, shared_yaxes=True)
    fig.add_trace(fig1.data[0], row=1, col=1)
    fig.add_trace(fig2.data[0], row=1, col=2)
    fig.update_layout(coloraxis={'showscale': False, 'colorscale': 'Inferno'}, showlegend=False, xaxis={'autorange': "reversed"},
                      yaxis={'autorange': True}, xaxis2={'autorange': "reversed"}, autosize=False)
    fig.update_xaxes(title_text="RA (mas)", constrain="domain")
    fig.update_yaxes(title_text="Dec (mas)", row=1, col=1, scaleanchor="x", scaleratio=1)

    return fig
    return fig1
            # layout :  xaxis autorange reversed





#####################  This is the webpage layout
app.layout = html.Div([
    html.Div(id='banner', className='navbar-brand d-flex p-3 shadow-sm', children=[
        html.A(className='d-inline-block mr-md-auto', href="https://www.evlbi.org", children=[
            html.Img(height='70px', src=app.get_asset_url("logo_evn.png"),
                     alt='European VLBI Network (EVN)',
                     className="d-inline-block align-top"),
        ]),
        html.H2('EVN Observation Planner', className='d-inline-block align-middle mx-auto'),
        html.A(className='d-inline-block ml-auto pull-right', href="https://www.jive.eu", children=[
            html.Img(src=app.get_asset_url("logo_jive.png"), height='70px',
                     alt='Joinst Institute for VLBI ERIC (JIVE)')
        ])
    ]),
    html.Div([html.Br(), html.Br()]),
    # html.Div(id='full-window', children=html.Div(id='main-window', children=main_page(False)))])
    html.Div(id='full-window', children=initial_page())])





if __name__ == '__main__':
    # app.run_server(host='0.0.0.0', debug=True)
    # app.run_server(debug=True)
    app.run_server(host='0.0.0.0', debug=True)



