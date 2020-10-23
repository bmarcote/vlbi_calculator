#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""EVN Observation Planner.

Program to compute the source elevation visibility
and expected thermal noise level for a given EVN observation.
"""

__author__ = "Benito Marcote"
__credits__ = "Benito Marcote"
__license__ = "LGPLv3+"
__date__ = "2020/04/21"
__version__ = "1.0rc1"
__maintainer__ = "Benito Marcote"
__email__ = "marcote@jive.eu"
__status__ = "Development"   # Prototype, Development, Production.

import os
from os import path
from time import sleep
import itertools
from importlib import resources
from datetime import datetime as dt
import numpy as np
import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import plotly.graph_objs as go
from astropy.time import Time
from astropy import coordinates as coord
from astropy import units as u

## THIS WILL NEED TO GO AWAY IN THE NEW VERSION OF ASTROPY, WHICH IS STILL NOT
## SUPPORTED BY THE CURRENT VERSION OF ASTROPLAN
# Tweak to not let astroplan crashing...
from astropy.utils.data import clear_download_cache
from astropy.utils import iers
clear_download_cache()  # to be sure it is really working

iers.conf.auto_download = False
iers.conf.iers_auto_url = None
iers.conf.auto_max_age = None
iers.conf.remote_timeout = 100.0
iers.conf.download_cache_lock_attempts = 10

from astroplan import FixedTarget

from vlbiplanobs import freqsetups as fs
from vlbiplanobs import stations
from vlbiplanobs import observation
from vlbiplanobs import graphical_elements as ge
# adding the possibility of disabled. Will be implemented in a future version of dash_bootstrap_components
from vlbiplanobs.Checkbox import Checkbox


all_antennas = stations.Stations.get_stations_from_configfile()

sorted_networks = {'EVN': 'EVN: European VLBI Network', 'eMERLIN': 'eMERLIN (out-stations)',
                   'VLBA': 'VLBA: Very Long Baseline Array',
                   'LBA': 'LBA: Australian Long Baseline Array',
                   'KVN': 'KVN: Korean VLBI Network',
                   'Other': 'Other antennas'}
default_arrays = {'EVN': ['Ef', 'Hh', 'Jb2', 'Mc', 'Nt', 'Ur', 'On', 'Sr', 'T6', 'Tr',
                          'Ys', 'Wb', 'Bd', 'Sv', 'Zc', 'Ir'],
          'e-EVN': ['Ef', 'Hh', 'Jb2', 'Mc', 'Nt', 'On', 'T6', 'Tr', 'Ys', 'Wb',
                    'Bd', 'Sv', 'Zc', 'Ir', 'Sr', 'Ur'],
          'eMERLIN': ['Cm', 'Kn', 'Pi', 'Da', 'De', 'Jb2'],
          'LBA': ['ATCA', 'Pa', 'Mo', 'Ho', 'Cd', 'Td', 'Ww'],
          'VLBA': ['Br', 'Fd', 'Hn', 'Kp', 'La', 'Mk', 'Nl', 'Ov', 'Pt', 'Sc'],
          'KVN': ['Ky', 'Ku', 'Kt'],
          'Global VLBI': ['Ef', 'Hh', 'Jb2', 'Mc', 'Nt', 'Ur', 'On', 'Sr', 'T6',
                          'Tr', 'Ys', 'Wb', 'Bd', 'Sv', 'Zc', 'Ir', 'Br', 'Fd', 'Hn',
                          'Kp', 'La', 'Mk', 'Nl', 'Ov', 'Pt', 'Sc'],
          'HSA': ['Br', 'Fd', 'Hn', 'Kp', 'La', 'Mk', 'Nl', 'Ov', 'Pt', 'Sc', 'Ef', 'Ar',
                  'Gb', 'Y27'],
          'GMVA': ['Ef', 'Mh', 'On', 'Ys', 'Pv', 'Br', 'Fd', 'Kp', 'La', 'Mk', 'Nl',
                   'Ov', 'Pt'],
          'EHT': ['ALMA', 'Pv', 'LMT', 'PdB', 'SMA', 'JCMT', 'APEX', 'SMT', 'SPT']}
default_datarates = {'EVN': 2048, 'e-EVN': 2048, 'eMERLIN': 4096, 'LBA': 1024, 'VLBA': 4096, 'KVN': 4096,
                     'Global VLBI': 2048, 'HSA': 2048, 'GMVA': 4096, 'EHT': 2**15}


# Safety check that all these antennas are available in the file
for a_array in default_arrays:
    for a_station in default_arrays[a_array]:
        assert a_station in all_antennas.keys()

doc_files = {'About this tool': 'doc-contact.md',
             'About the antennas': 'doc-antennas.md',
             'Technical background': 'doc-estimations.md'}

# Initial values
# target_source = observation.Source('10h2m3s +50d40m30s', 'Source')
# obs_times = Time('1967-04-17 10:00') + np.arange(0, 600, 15)*u.min

selected_band = ''
obs = observation.Observation()
# obs = observation.Observation(target=target_source)
# obs.times = Time('2020-06-15 20:00', scale='utc') + np.arange(0, 1200, 30)*u.min
# obs.band = selected_band
# obs.datarate = 1024
# obs.subbands = 8
# obs.channels = 32
# obs.polarizations = 2
# obs.inttime = 2




# external_stylesheets = ["https://stackpath.bootstrapcdn.com/bootstrap/4.5.0/css/bootstrap.min.css", "http://jive.eu/~marcote/style.css"]
# external_stylesheets = ["https://bmarcote.github.io/temp/style.css"]
external_stylesheets = []
# n_timestamps = 70 # Number of points (timestamps) for the whole observations.
external_scripts = ["https://kit.fontawesome.com/69c65a0ab5.js"]
#external_scripts = ["""https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js", \
#         "https://polyfill.io/v3/polyfill.min.js?features=es6"]
#         "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"]


app = dash.Dash(__name__, external_scripts=external_scripts)
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
                   [dcc.Markdown(parsed_text), ge.antenna_cards(app, all_antennas)], id=str(i))]
            else:
                temp += [ge.create_accordion_card(a_topic, dcc.Markdown(parsed_text),
                         id=str(i))]

    return html.Div(temp, className='col-12 accordion')


@app.callback([Output(f"collapse-{i}", "is_open") for i in range(len(doc_files))],
        [Input(f"group-{i}-toggle", "n_clicks") for i in range(len(doc_files))],
        [State(f"collapse-{i}", "is_open") for i in range(len(doc_files))])
def toggle_accordion(*args):
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
    """Message written in a modal error window.
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

    return colon_coord.replace(' ', 's ')+'s'


def alert_message(message, title="Warning!"):
    """Produces an alert-warning message.
    message can be either a string or a list with different string/dash components.
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
    """Given the observation, it sets the text about all properties of the observation.
    """
    cards = []
    # The time card
    cards += ge.summary_card_times(app, obs)
    cards += ge.summary_card_antennas(app, obs)
    cards += ge.summary_card_beam(app, obs)
    cards += ge.summary_card_frequency(app, obs)
    cards += ge.summary_card_rms(app, obs)
    cards += ge.summary_card_fov(app, obs)
    return [html.Div(className='card-deck col-12 justify-content-center', children=cards)]



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
    html.Div([html.Br()]),
    # First row containing all buttons/options, list of telescopes, and button with text output
    dcc.ConfirmDialog(id='global-error', message=''),
    # Elements in second column (checkboxes with all stations)
    html.Div(className='container-fluid', children=[
        dcc.Tabs(parent_className='custom-tabs', className='custom-tabs-container', children=[
            dcc.Tab(label='Observation Setup', className='custom-tab',
                    selected_className='custom-tab--selected', children=[
                # Elements in first column ()
                html.Div(className='row justify-content-center', children=[
                html.Div(className='col-sm-3', style={'max-width': '350px','float': 'left',
                                                      'padding': '2%'}, children=[
                    html.Div(className='form-group', children=[
                        html.Label('Select default VLBI Network(s)'),
                        *ge.tooltip(idname='popover-network', message="Automatically selects "
                                    "the default antennas for the selected VLBI network(s)."),
                        dcc.Dropdown(id='array', options=[{'label': n, 'value': n} \
                                for n in default_arrays if n != 'e-EVN'], value=[],
                                multi=True),
                    ]),
                    html.Div(className='input-group-prepend', children=[
                        dbc.Checklist(id='e-EVN', className='checkbox', persistence=True,
                                      options=[{'label': ' e-EVN (real-time) mode',
                                                'value': 'e-EVN'}], value=[]),
                        *ge.tooltip(idname='popover-eevn',
                            message="Only available for the EVN: real-time correlation mode.")
                    ]),
                    html.Br(),
                    html.Div(className='form-group', children=[
                        html.Label('Start of observation (UTC)'),
                        *ge.tooltip(idname='popover-startime', message="Select the date and "
                                    "time of the start of the observation (Universal, UTC, "
                                    "time). You will also see the day of the year (DOY) in "
                                    "brackets once the date is selected."),
                        dcc.DatePickerSingle(id='starttime', min_date_allowed=dt(1900, 1, 1),
                                             max_date_allowed=dt(2100, 1, 1),
                                             display_format='DD-MM-YYYY (DDD)',
                                             placeholder='Start date',
                                             first_day_of_week=1,
                                             initial_visible_month=dt.today(),
                                             persistence=True,
                                             className='form-picker'),
                        dcc.Dropdown(id='starthour', placeholder="Start time",
                                     options=[{'label': f"{hm//60:02n}:{hm % 60:02n}", \
                                               'value': f"{hm//60:02n}:{hm % 60:02n}"} \
                                              for hm in range(0, 24*60, 15)],
                                     persistence=True, className='form-hour'),
                        html.Small(id='error_starttime', style={'color': 'red'},
                                   className='form-text text-muted')
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('End of observation (UTC)'),
                        *ge.tooltip(idname='popover-endtime', message="Select the date and "
                                    "time of the end of the observation (Universal, UTC, "
                                    "time). You will also see the day of the year (DOY) in "
                                    "brackets once the date is selected."),
                        dcc.DatePickerSingle(id='endtime', min_date_allowed=dt(1900, 1, 1),
                                             max_date_allowed=dt(2100, 1, 1),
                                             display_format='DD-MM-YYYY (DDD)',
                                             placeholder='End date',
                                             first_day_of_week=1,
                                             initial_visible_month=dt.today(),
                                             persistence=True,
                                             className='form-picker'),
                        dcc.Dropdown(id='endhour', placeholder="End time",
                                     options=[{'label': f"{hm//60:02n}:{hm % 60:02n}", \
                                               'value': f"{hm//60:02n}:{hm % 60:02n}"} \
                                              for hm in range(0, 24*60, 15)],
                                     persistence=True, className='form-hour'),
                        html.Small(id='error_endtime', style={'color': 'red'},
                                   className='form-text text-muted')
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('Target Source Coordinates'),
                        *ge.tooltip(idname='popover-target',
                                 message="Right Ascension and Declination of the source " \
                                         "J2000 coordinates are assumed. Both, 00:00:00 00:00:00 and " \
                                         "00h00m00s 00d00m00s syntaxes are allowed."),
                        # dcc.Input(id='source', value='12:29:06.7 +02:03:08.6', type='text',
                        dcc.Input(id='source', value=None, type='text',
                                  className='form-control', placeholder="hh:mm:ss dd:mm:ss",
                                  persistence=True),
                        html.Small(id='error_source', style={'color': 'red'},
                                   className='form-text text-muted'),
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label(id='onsourcetime-label',
                                   children='% of on-target time'),
                        *ge.tooltip(idname='popover-ontarget',
                                 message="Assumes that you will only spend this amount of the total " \
                                         "observing time on the given target source. It affects the " \
                                         "expected sensitivity."),
                        dcc.Slider(id='onsourcetime', min=20, max=100, step=5, value=75,
                                   marks= {i: str(i) for i in range(20, 101, 10)},
                                   persistence=True),
                    ]),
                    html.H4("Advanced setup"),
                    html.Div(className='form-group', children=[
                        html.Label('Datarate per station'),
                        *ge.tooltip(idname='popover-datarate',
                                 message=["Expected datarate for each station, assuming all " \
                                          "of them run at the same rate.",
                                     html.Ul([
                                        html.Li("The EVN can run typically at up to 2 Gbps (1 Gbps at L band), " \
                                                "although a few antennas may observe at lower datarates."),
                                        html.Li("The VLBA can now observe up to 4 Gbps."),
                                        html.Li("The LBA typically runs at 512 Mbps but can reach up to 1 Gbps."),
                                        html.Li("Check the documentation from other networks to be " \
                                                "sure about their capabilities.")])]),
                        dcc.Dropdown(id='datarate',
                                     placeholder="Select the data rate...",
                                     options=[{'label': fs.data_rates[dr], 'value': dr} \
                                     for dr in fs.data_rates], value=2048, persistence=True),
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('Number of subbands'),
                        *ge.tooltip(idname='popover-subbands',
                                 message="Number of subbands the total observed bandwidth "
                                         "will be split during correlation."),
                        dcc.Dropdown(id='subbands', placeholder="Select no. subbands...",
                                     options=[{'label': f"{sb} subbands", 'value': sb} \
                                     for sb in fs.subbands], value=8, persistence=True),
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('Number of spectral channels'),
                        *ge.tooltip(idname='popover-channels',
                                 message="How many channels per subband will be produced "
                                         "during correlation."),
                            dcc.Dropdown(id='channels', placeholder="Select no. channels...",
                                         options=[{'label': f"{ch} channels per subband",
                                                   'value': ch} \
                                         for ch in fs.channels], value=32, persistence=True),
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('Number of polarizations'),
                        *ge.tooltip(idname='popover-pols',
                            message="Number of polarizations to correlate. Note that VLBI uses circular " \
                                    "polarizations. Full polarization implies the four stokes: RR, LL, RL, LR; " \
                                    "while dual polarization implies RR and LL only."),
                        dcc.Dropdown(id='pols', placeholder="Select polarizations...",
                                     options=[{'label': fs.polarizations[p], 'value': p} \
                                     for p in fs.polarizations], value=4, persistence=True),
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('Integration time (s)'),
                        *ge.tooltip(idname='popover-inttime',
                            message="Integration time to compute each visibility. Note that for continuum " \
                                    "observations values of 1-2 seconds are typical."),
                        dcc.Dropdown(id='inttime', placeholder="Select integration time...",
                                     options=[{'label': fs.inttimes[it], 'value': it} \
                                     for it in fs.inttimes], value=2, persistence=True),
                    ])
                ]),
                html.Div(className='col-9', children=[
                    html.Div(id='first-advise', className='col-sm-9', children=[
                        html.P(["This is the ", html.B("EVN Observation Planner"),
                               ". First select the "
                               "band (frequency/wavelength) at which you want to observe, "
                               "then customize your observation setup (left options and "
                               "select wished antennas), and finally "
                               "run ", html.B("'compute observation'"),
                               ". You will get a detailed "
                               "summary of the planned observation (like when the source "
                               "is visible, expected rms noise level, etc.) in the different "
                               "tabs."]),
                    ], style={'margin-top': '2rem', 'margin-bottom': '2rem'}),
                    html.Div(className='col-9 form-group row align-items-end', children=[
                        html.Div(className='col-md-6', children=[
                            html.Label('First Step: Select your observing Band',
                                        style={'color': '#a01d26'}),
                            *ge.tooltip(idname='popover-band',
                                    message="First select at which frequency/wavelength "
                                            "you want to observe. This will update the "
                                            "antenna list showing the ones that can observe "
                                            "at that given frequency."),
                            dcc.Dropdown(id='band', persistence=True,
                                 options=[{'label': fs.bands[b], 'value': b} for b \
                                # in fs.bands], value='18cm'),
                                in fs.bands], placeholder='Select observing band...')
                        ]),
                        html.Div(className='col-sm-3', children=[
                            html.Button('Compute Observation', id='antenna-selection-button',
                                        className='btn btn-primary btn-lg'),
                        ])
                    ]),
                    html.Div(className='col-9 text-center justify-content-center', children=[
                        dcc.Loading(id="loading", children=[html.Div(id="loading-output")],
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
                                             className='custom-control-label form-check-label')
                                ], check=True, inline=True, className="form-check-input "
                                "custom-checkbox custom-control custom-control-inline")
                                for s in all_antennas if s.network == an_array
                            ])
                        ]) for an_array in sorted_networks
                    ])
                ]),
                # html.Div(className='col-sm-2', style={'float': 'left'}, children=[
                # ])
                ])
            ]),

            dcc.Tab(label='Summary', className='custom-tab',
                    selected_className='custom-tab--selected', children=[
                html.Div(className='row justify-content-center', children=[
                    html.Div(className='col-10 justify-content-center',
                             id='sensitivity-output',
                             children=[html.Div(className='col-md-6', children=[
                                html.Br(), html.Br(), html.H2("Set the observation first"),
                                html.P("Here you will see a summary of your observation, "
                                       "with information about all participating stations, longest and "
                                       "shortest baseline, expected size of the data once is correlated, "
                                       "reached resolution and sensitivity, and the limitations in your "
                                       "field of view due to time and frequency smearing.")])
                            ])
                ])
            ]),
            dcc.Tab(label='Elevations', className='custom-tab',
                    selected_className='custom-tab--selected', children=[
                html.Div(className='row justify-content-center', children=[
                html.Div(className='col-md-8 justify-content-center', children=[
                    # Elevation VS time
                    html.Br(),
                    html.Div([
                        html.Br(),
                        html.H4("When is your source visible?"),
                        html.Br(),
                        dbc.Alert([html.H4("Info on plots", className='alert-heading'),
                                   html.P("A single click on one station in the legend will "
                                           "hide/show it. Double-click will hide/show "
                                           "all other antennas. You can also save the plot "
                                           "as png."),
                                  ], color='info', dismissable=True),
                        html.Br(),
                        html.P("The following plot shows the source elevation for the "
                        "different antennas during the proposed observation. The horizontal "
                        "solid and dashed lines represent the elevation of 20 and 10 degrees, "
                        "respectively.")
                    ]),
                    html.Div([
                        dcc.Graph(id='fig-elev-time')
                    ],className='tex2jax_ignore'),
                    html.Br(),
                    html.Div([
                        html.P("""The following plot shows when the source may be observed
                        for the different antennas, assuming a minimum elevation of 10 degrees
                        for most antennas (except e.g. Arecibo). Note that some antennas may
                        have additional constraints por some particular azimuth or elevation
                        angles that are not considered here.
                        """)
                    ]),
                    html.Div([
                        dcc.Graph(id='fig-ant-time')
                    ],className='tex2jax_ignore')
                ])])
            ]),
            dcc.Tab(label='Coverage', className='custom-tab',
                    selected_className='custom-tab--selected', children=[
                #  Images
                html.Div(className='row justify-content-center', children=[
                html.Div(className='col-md-8 justify-content-center', children=[
                    # dcc.Markdown(children="""To be implemented.
                    #     The uv coverage and expected dirty images will go here.""")
                    html.Div(className='col-md-8 justify-content-center',
                             children=[dcc.Graph(id='fig-uvplane')])
                ])])
            ]),
            dcc.Tab(label='Documentation', className='custom-tab',
                    selected_className='custom-tab--selected', children=[
                #  Documentation
                html.Div(className='row justify-content-center', children=[
                    html.Div([html.Br(), html.Br()]),
                    html.Div(className='col-md-8', children=get_doc_text())
                ])
            ])
        ]),
    html.Div(className='container-fluid', children=[html.Br(), html.Br()])
    ])
])





@app.callback(Output('onsourcetime-label', 'children'),
              [Input('onsourcetime', 'value')])
def update_onsourcetime_label(onsourcetime):
    return f"% of on-target time ({onsourcetime}%)"


@app.callback([Output(f"check_{s.codename}", 'checked') for s in all_antennas] + \
              [Output(f"check_{s.codename}", 'disabled') for s in all_antennas] + \
              [Output('datarate', 'value')],
              [Input('band', 'value'), Input('array', 'value'), Input('e-EVN', 'value')])
def select_antennas(selected_band, selected_networks, is_eEVN):
    """Given a selected band and selected default networks, it selects the associated
    antennas from the antenna list.
    """
    selected_antennas = []
    if is_eEVN:
        selected_antennas = [ant for ant in default_arrays['e-EVN'] \
                             if all_antennas[ant].has_band(selected_band)]
        datarate = default_datarates['e-EVN'] if selected_band not in ('18cm', '21cm') else 1024
        return [True if s.codename in selected_antennas else False for s in all_antennas] + \
               [False if (s.has_band(selected_band) and s.real_time) else True \
                for s in all_antennas] + [datarate]
    else:
        datarate = -1
        for an_array in selected_networks:
            selected_antennas += [ant for ant in default_arrays[an_array] \
                                    if all_antennas[ant].has_band(selected_band)]
            datarate = max(datarate, default_datarates[an_array] if not ((an_array == 'EVN') and \
                                                    (selected_band in ('18cm', '21cm'))) else 1024)

        return [True if s.codename in selected_antennas else False for s in all_antennas] + \
               [False if s.has_band(selected_band) else True for s in all_antennas] + [datarate]



@app.callback([Output('error_starttime', 'children'),
               Output('error_endtime', 'children')],
              [Input('starttime', 'date'), Input('starthour', 'value'),
               Input('endtime', 'date'), Input('endhour', 'value')])
def check_obstime(starttime, starthour, endtime, endhour):
    times = [None, None]
    if None not in (starttime, starthour):
        times[0] = Time(dt.strptime(f"{starttime} {starthour}", '%Y-%m-%d %H:%M'))

    if None not in (endtime, endhour):
        times[1] = Time(dt.strptime(f"{endtime} {endhour}", '%Y-%m-%d %H:%M'))

    if None in times:
        return '', ''

    if (times[1] - times[0]) > 5*u.d:
        return ["Please, put an observation shorter than 5 d"]*2
    elif (times[0] - times[1]) >= 0*u.d:
        return ["Start time must be earlier than end time"]*2

    return '', ''


@app.callback(Output('error_source', 'children'),
        [Input('source', 'value')])
def get_source(source_coord):
    if source_coord != 'hh:mm:ss dd:mm:ss':
        try:
            dummy_target = observation.Source(convert_colon_coord(source_coord), 'Source')
            return ''
        except ValueError as e:
            return "Use 'hh:mm:ss dd:mm:ss' format"
    else:
        return dash.no_update



@app.callback([Output('loading-output', 'children'),
               Output('sensitivity-output', 'children'),
               Output('fig-elev-time', 'figure'),
               Output('fig-ant-time', 'figure'),
               Output('fig-uvplane', 'figure'), Output('global-error', 'message')],
              [Input('antenna-selection-button', 'n_clicks')],
              [State('band', 'value'),
               State('starttime', 'date'),
               State('starthour', 'value'),
               State('endtime', 'date'),
               State('endhour', 'value'),
               State('source', 'value'),
               State('onsourcetime', 'value'),
               State('datarate', 'value'),
               State('subbands', 'value'),
               State('channels', 'value'),
               State('pols', 'value'),
               State('inttime', 'value')] + \
               [State(f"check_{s.codename}", 'checked') for s in all_antennas])
def compute_observation(n_clicks, band, starttime, starthour, endtime, endhour, source,
                        onsourcetime, datarate, subbands, channels, pols, inttime, *ants):
    """Computes all products to be shown concerning the set observation.
    """
    if n_clicks is None:
        return '', dash.no_update, dash.no_update, dash.no_update, \
               dash.no_update, dash.no_update
    try:
        target_source = observation.Source(convert_colon_coord(source), 'Source')
    except ValueError as e:
        return alert_message(["Incorrect format for source coordinates.", html.Br(),
                f"{'Empty value' if source=='' else source} found but 'hh:mm:ss dd:mm:ss' expected."]),\
               "First, set correctly an observation in the previous tab.", \
               dash.no_update, dash.no_update, dash.no_update, dash.no_update
    try:
        time0 = Time(dt.strptime(f"{starttime} {starthour}", '%Y-%m-%d %H:%M'),
                     format='datetime', scale='utc')
    except ValueError as e:
        return alert_message("Incorrect format for starttime."), \
               "First, set correctly an observation in the previous tab.", \
               dash.no_update, dash.no_update, dash.no_update, dash.no_update

    try:
        time1 = Time(dt.strptime(f"{endtime} {endhour}", '%Y-%m-%d %H:%M'),
                     format='datetime', scale='utc')
    except ValueError as e:
        return alert_message("Incorrect format for endtime."), \
               "First, set correctly an observation in the previous tab.", \
               dash.no_update, dash.no_update, dash.no_update, dash.no_update

    if time0 >= time1:
        return alert_message("The start time must be earlier than the end time"), \
               "First, set correctly an observation in the previous tab.", \
               dash.no_update, dash.no_update, dash.no_update, dash.no_update

    if (time1 - time0) > 5*u.d:
        return alert_message("Please, set an observation that last for less than 5 days."), \
               "First, set correctly an observation in the previous tab.", \
               dash.no_update, dash.no_update, dash.no_update, dash.no_update

    # try:
    # TODO: this should not be hardcoded...
    obs_times = time0 + np.linspace(0, (time1-time0).to(u.min).value, 50)*u.min
    # all_selected_antennas = list(itertools.chain.from_iterable(ants))
    try:
        obs = observation.Observation(target=target_source, times=obs_times, band=band,
                      datarate=datarate, subbands=subbands, channels=channels,
                      polarizations=pols, inttime=inttime, ontarget=onsourcetime/100.0,
                      stations=stations.Stations('Observation',
                                                 itertools.compress(all_antennas, ants)))
        sensitivity_results = update_sensitivity(obs)
    except observation.SourceNotVisible:
        return alert_message([
                    html.P(["Your source cannot be observed within the arranged observation.",
                    html.Br(),
                    "There are no antennas that can simultaneously observe your source "
                    "during the given observing time."]),
                    html.P("Modify the observing time or change the selected antennas"
                           " to observe this source.")], title="Warning!"), \
                dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

    # return update_sensitivity(obs), dash.no_update, dash.no_update
    # TODO: parallelize all these functions
    return [html.Br(),
            dbc.Alert("You can check now the results in the different tabs", color='info', \
                      dismissable=True)], sensitivity_results, \
           get_fig_ant_elev(obs), get_fig_ant_up(obs), get_fig_uvplane(obs), dash.no_update



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
    data_fig.append({'x': np.unwrap(obs.gstimes.value), 'y': np.zeros_like(obs.times)+20,
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
                                            np.floor(np.unwrap(obs.gstimes.value)[-1])+1),
                                  'ticktext': np.arange(np.ceil(obs.gstimes.value[0]),
                                            np.floor(np.unwrap(obs.gstimes.value)[-1])+1) % 24,
                                  'ticks': 'inside', 'showline': True, 'mirror': False,
                                  'hovermode': 'closest', 'color': 'black', 'side': 'top'},
                       'yaxis': {'title': 'Elevation (degrees)', 'range': [0., 92.],
                                 'ticks': 'inside', 'showline': True, 'mirror': "all",
                                 'showgrid': False, 'hovermode': 'closest'},
                                 'zeroline': True, 'zerolinecolor': 'k'}}



def get_fig_ant_up(obs):
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
            'layout': {'title': 'uv coverage', 'showlegend': False,
                       'hovermode': 'closest',
                       'width': 700, 'height': 700,
                       'xaxis': {'title': 'u (lambda)', 'showgrid': False, 'zeroline': False,
                                 'ticks': 'inside', 'showline': True, 'mirror': "all",
                                 'color': 'black'},
                       'yaxis': {'title': 'v (lambda)', 'showgrid': False, 'scaleanchor': 'x',
                                 'ticks': 'inside', 'showline': True, 'mirror': "all",
                                 'color': 'black', 'zeroline': False}}}


def get_fig_dirty_map(obs):
    pass




if __name__ == '__main__':
    # app.run_server(host='0.0.0.0', debug=True)
    # app.run_server(debug=True)
    app.run_server(host='0.0.0.0')


