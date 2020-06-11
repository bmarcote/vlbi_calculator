#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""EVN Calculator.

Program to compute the source elevation visibility and expected thermal
noise level for a given EVN observation.
"""

__author__ = "Benito Marcote"
__copyright__ = "Copyright 2020, Joint Insitute for VLBI-ERIC (JIVE)"
__credits__ = "Benito Marcote"
__license__ = "GPL"
__date__ = "2020/04/21"
__version__ = "0.0.1"
__maintainer__ = "Benito Marcote"
__email__ = "marcote@jive.eu"
__status__ = "Development"   # Prototype, Development, Production.

from os import path
from time import sleep
import itertools
import datetime
import numpy as np
import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import plotly.graph_objs as go
from datetime import datetime as dt
from astropy.time import Time
from astropy import coordinates as coord
from astropy import units as u
# Tweak to not let astroplan crashing...

from astropy.utils.data import clear_download_cache
clear_download_cache()  # to be sure it is really working
from astropy.utils import iers
iers.conf.auto_download = False
iers.conf.iers_auto_url = None

from astroplan import FixedTarget
from src import freqsetups as fs
from src import stations
from src import functions as fx
from src import observation


current_directory = path.dirname(path.realpath(__file__))
# stationList =  stations.Stations()
# stationList.add_from_file(current_directory+'/station_location.txt')


# iers.IERS.iers_table = iers.IERS.open(cache=True)
# iers.IERS.iers_table = iers.IERS_A.open(iers.IERS_A_URL)
# iers.Conf.iers_auto_url.set('ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all')

all_antennas = fx.get_stations_from_file(f"{current_directory}/data/station_location.txt")
sorted_networks = {'EVN': 'EVN: European VLBI Network', 'eMERLIN': 'eMERLIN',
                   'VLBA': 'VLBA: Very Long Baseline Array',
                   'LBA': 'LBA: Australian Long Baseline Array',
                   'KVN': 'Korean VLBI Network',
                   'Other': 'Other antennas'}
default_arrays = {'EVN': ['Ef', 'Hh', 'Jb2', 'Mc', 'Nt', 'Ur', 'On', 'Sr', 'T6', 'Tr',
                          'Ys', 'Wb', 'Bd', 'Sv', 'Zc', 'Ir'],
          'e-EVN': ['Ef', 'Hh', 'Ir', 'Jb2', 'Mc', 'Nt', 'On', 'T6', 'Tr', 'Ys', 'Wb',
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


# Safety check that all these antennas are available in the file
for a_array in default_arrays:
    for a_station in default_arrays[a_array]:
        assert a_station in all_antennas.keys()


# Initial values
target_source = observation.Source('1h2m3s +50d40m30s', 'Source')
# obs_times = Time('1967-04-17 10:00') + np.arange(0, 600, 15)*u.min
selected_band = '18cm'




# external_stylesheets = ["https://stackpath.bootstrapcdn.com/bootstrap/4.5.0/css/bootstrap.min.css", "http://jive.eu/~marcote/style.css"]
# external_stylesheets = ["https://bmarcote.github.io/temp/style.css"]
external_stylesheets = []
# n_timestamps = 70 # Number of points (timestamps) for the whole observations.





# app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app = dash.Dash(__name__)


server = app.server

# app.config.requests_pathname_prefix = ''
# app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})

# app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/brPBPO.css"})

def tooltip(message, idname, trigger='?', placement='right', **kwargs):
    """Defines a tooltip (popover) that will be shown in the rendered page.
    It will place a <sup>`trigger`</sup> in the page within a span tag.
    Returns the list with all html elements.
    """
    return [html.Span(children=html.Sup(trigger, className='popover-link'), id=idname),
            dbc.Tooltip(message, target=idname, placement=placement,
                        # className='tooltip-class', #innerClassName='tooltip-class-inner',
                        **kwargs)]



#####################  This is the webpage layout
app.layout = html.Div([
    html.Script("""
 $(document).ready(function() {
  $('[data-toggle="popover"]').popover({
      placement : 'bottom',
      title : '<div style="text-align:center; color:red; text-decoration:underline; font-size:14px;"> Muah ha ha</div>', //this is the top title bar of the popover. add some basic css
      html: 'true',
      content : '<div id="popOverBox"><img src="http://www.hd-report.com/wp-content/uploads/2008/08/mr-evil.jpg" width="251" height="201" /></div>' //this is the content of the html box. add the image here or anything you want really.
});
});
"""),
    html.Div(id='banner', className='navbar-brand d-flex p-3 shadow-sm', children=[
        html.A(className='d-inline-block mr-md-auto', href="https://www.evlbi.org", children=[
            html.Img(height='70px', src=app.get_asset_url("logo_evn.png"),
                     alt='European VLBI Network (EVN)',
                     className="d-inline-block align-top"),
        ]),
        html.H2('EVN Source Visibility', className='d-inline-block align-middle mx-auto'),
        html.A(className='d-inline-block ml-auto pull-right', href="https://www.jive.eu", children=[
            html.Img(src=app.get_asset_url("logo_jive.png"), height='70px',
                     alt='Joinst Institute for VLBI ERIC (JIVE)')
        ])
        # html.Img(src='http://www.ira.inaf.it/evnnews/archive/evn.gif')
    ]),
    # ], className='banner'),
    html.Div([html.Br()]), #style={'clear': 'both', 'margin-top': '20px'}),
    # First row containing all buttons/options, list of telescopes, and button with text output
    dcc.ConfirmDialog(id='global-error', message=''),
    # Elements in second column (checkboxes with all stations)
    html.Div(className='container-fluid', children=[
        dcc.Tabs(parent_className='custom-tabs', className='custom-tabs-container', children=[
            dcc.Tab(label='Observation Setup', className='custom-tab',
                    selected_className='custom-tab--selected', children=[
                # Elements in first column ()
                html.Div(className='row justify-content-center', children=[
                html.Div(className='col-sm-3', style={'max-width': '300px','float': 'left',
                                                      'padding': '2%'}, children=[
                    html.Div(className='form-group', children=[
                        html.Label('Observing Band'),
                        *tooltip(idname='popover-band', message="First select the observing band. Antenna list will be updated and only the ones that can observe at this band will be enable."),
                        dcc.Dropdown(id='band', persistence=True,
                                 options=[{'label': fs.bands[b], 'value': b} for b \
                                in fs.bands], value='18cm'),
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('Select default VLBI Network(s)'),
                        *tooltip(idname='popover-network', message="Automatically selects the default antennas for the selected VLBI network(s)."),
                        dcc.Dropdown(id='array', options=[{'label': n, 'value': n} \
                                for n in default_arrays if n != 'e-EVN'], value=['EVN'],
                                multi=True),
                    ]),
                    html.Div(className='input-group-prepend', children=[
                        dcc.Checklist(id='e-EVN', className='checkbox', persistence=True,
                                      options=[{'label': ' e-EVN (real-time) mode',
                                                'value': 'e-EVN'}], value=[]),
                        *tooltip(idname='popover-eevn',
                            message="Only available for the EVN: real-time correlation mode.")
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('Start of observation (UTC)'),
                        # dcc.Input(id='starttime', value='DD/MM/YYYY HH:MM', type='text',
                        dcc.Input(id='starttime', value='17/04/1967 10:00', type='text',
                                  className='form-control', placeholder="dd/mm/yyyy HH:MM",
                                  persistence=True),
                        html.Small(id='error_starttime', style={'color': 'red'},
                                   className='form-text text-muted')
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('End of observation (UTC)'),
                        # dcc.Input(id='endtime', value='DD/MM/YYYY HH:MM', type='text',
                        dcc.Input(id='endtime', value='17/04/1967 20:00', type='text',
                                  className='form-control', placeholder="dd/mm/yyyy HH:MM",
                                  persistence=True),
                        html.Small(id='error_endtime', style={'color': 'red'},
                                   className='form-text text-muted')
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('Target Source Coordinates'),
                        *tooltip(idname='popover-target',
                                 message="J2000 coordinates are assumed."),
                        # dcc.Input(id='source', value='hh:mm:ss dd:mm:ss', type='text',
                        dcc.Input(id='source', value='12:29:06.7 +02:03:08.6', type='text',
                                  className='form-control', placeholder="hh:mm:ss dd:mm:ss",
                                  persistence=True),
                        html.Small(id='error_source', style={'color': 'red'},
                                   className='form-text text-muted'),
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label(id='onsourcetime-label',
                                   children='% of on-target time'),
                        *tooltip(idname='popover-ontarget',
                                 message="Assumes that you will only spend this amount of the total observing time on the given target source. It affects the expected sensitivity."),
                        dcc.Slider(id='onsourcetime', min=20, max=100, step=5, value=75,
                                   marks= {i: str(i) for i in range(20, 101, 10)},
                                   persistence=True),
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('Datarate per station (in Mbps)'),
                        *tooltip(idname='popover-datarate',
                                 message=["Expected datarate for each station, assuming all of them run at the same rate.",
                                     html.Ul([
                                        html.Li("The EVN can run typically at up to 2 Gbps (1 Gbps at L band), although a few antennas may observe at lower datarates."),
                                        html.Li("The VLBA can now observe up to 4 Gbps."),
                                        html.Li("The LBA typically observes at 512 Mbps but can run up to 1 Gbps."),
                                        html.Li("Check the documentation from other networks to be sure about their capabilities.")])]),
                        dcc.Dropdown(id='datarate', placeholder="Select a datarate...",
                                     options=[{'label': str(dr), 'value': dr} \
                                     for dr in fs.data_rates], value=1024, persistence=True),
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('Number of subbands'),
                        *tooltip(idname='popover-subbands',
                                 message="In how many subbands the total band will be split during correlation."),
                        dcc.Dropdown(id='subbands', placeholder="Select no. subbands...",
                                     options=[{'label': str(sb), 'value': sb} \
                                     for sb in fs.subbands], value=8, persistence=True),
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('Number of spectral channels'),
                        *tooltip(idname='popover-channels',
                                 message="How many channels per subband will be produced after correlation."),
                            dcc.Dropdown(id='channels', placeholder="Select no. channels...",
                                         options=[{'label': str(ch), 'value': ch} \
                                         for ch in fs.channels], value=32, persistence=True),
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('Number of polarizations'),
                        *tooltip(idname='popover-pols',
                            message="Number of polarizations to correlate. Note that VLBI observes circular polarizations. Full polarization implies the four stokes: RR, LL, RL, LR; while dual polarization implies RR and LL only."),
                        dcc.Dropdown(id='pols', placeholder="Select polarizations...",
                                     options=[{'label': fs.polarizations[p], 'value': p} \
                                     for p in fs.polarizations], value=4, persistence=True),
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('Integration time (s)'),
                        *tooltip(idname='popover-inttime',
                            message="Integration time to compute each visibility. Note that for continuum observations values of 1-2 seconds are typical."),
                        dcc.Dropdown(id='inttime', placeholder="Select integration time...",
                                     options=[{'label': fs.inttimes[it], 'value': it} \
                                     for it in fs.inttimes], value=2, persistence=True),
                    ])
                ]),
                # html.Div(style={'margin-top': '20px'}, children=[
                html.Div(className='col-9', children=[
                    html.Div(className='col-9 text-center justify-content-center', children=[
                        html.Button('Compute Observation', id='antenna-selection-button',
                                className='btn btn-primary btn-lg'),
                        dcc.Loading(id="loading", children=[html.Div(id="loading-output")],
                                    type="dot")
                    ]),
                    html.Div(id='antennas-div', className='container', children=[
                    # List with all antennas
                        html.Div(className='antcheck', children=[html.Br(),
                            html.Label(html.H4(f"{sorted_networks[an_array]}")),
                            html.Br(),
                            dcc.Checklist(id=f"list_stations_{an_array}",
                                className='antcheck',
                                labelClassName='form-check-label',
                                inputClassName='form-check-input',
                                options=[{'label': s.name, 'value': s.codename,
                                'disabled': not s.has_band(selected_band)}
                                for s in all_antennas if s.network == an_array], value=[])
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
                # Sensitivity calculations
                    # dcc.Markdown(id='sensitivity-output',
                    #              children="Set the observation first.")
                    html.Div(className='col-12', id='sensitivity-output',
                            children=[html.Div(className='col-6 justify-content-center',
                                children=[html.Br(),
                                html.P("You need to set the observation and click in the 'Compute Observation' buttom first (go to the previous tab).")])
                    ])
                ])
            ]),
            dcc.Tab(label='Elevations', className='custom-tab',
                    selected_className='custom-tab--selected', children=[
                html.Div(className='row justify-content-center', children=[
                html.Div(className='col-md-8', children=[
                    # Elevation VS time
                    html.Div([
                        dcc.Graph(id='fig-elev-time')
                    ]),
                    # Antenna VS time (who can observe)
                    html.Div([
                        dcc.Graph(id='fig-ant-time')
                    ])
                ])])
            ]),
            dcc.Tab(label='Imaging', className='custom-tab',
                    selected_className='custom-tab--selected', children=[
                #  Images
                html.Div(className='row justify-content-center', children=[
                html.Div(className='col-md-8', children=[
                    # dcc.Markdown(children="""To be implemented.
                    #     The uv coverage and expected dirty images will go here.""")
                    html.Div([dcc.Graph(id='fig-uvplane')])
                ])])
            ]),
            dcc.Tab(label='Documentation', className='custom-tab',
                    selected_className='custom-tab--selected', children=[
                #  Documentation
                html.Div(className='row justify-content-center', children=[
                html.Div(className='col-md-8', children=[
                    dcc.Markdown(children="""The Help/doc.
                        All explanations and technical details will go here.""")
                ])])
            ])
        ])
    ])
])




def error_text(an_error):
    """Message written in a modal error window.
    """
    return f"An error occured.\n{an_error}.\nPlease report to marcote@jive.eu."



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


def get_selected_antennas(list_of_selected_antennas):
    """Given a list of antenna codenames, it returns a Stations object containing
    all given antennas.
    """
    selected_antennas = stations.Stations('Observation', [])

    for ant in list_of_selected_antennas:
        selected_antennas.add(all_antennas[ant])

    return selected_antennas


def optimal_units(value, units):
    """Given a value (with some units), returns the unit choice from all
    `units` possibilities that better suits the value.
    It is meant for the following use:
    Given 0.02*u.Jy and units = [u.kJy, u.Jy, u.mJy, u.uJy], it will
    return 20*u.mJy.
    units should have a decreasing-scale of units, and all of them
    compatible with the units in `value`.
    """
    for a_unit in units:
        if 0.8 < value.to(a_unit).value <= 800:
            return value.to(a_unit)

    # Value too high or too low
    if value.to(units[0]).value > 1:
        return value.to(units[0])

    return value.to(units[-1])


def create_sensitivity_card(title, message):
    """Defines one of the cards that are shown in the Sensitivity tab. Each tab
    shows a title `title` and a message `message`.
    If message is a list (of strings), it is assumed as different paragraphs.
    It returns the HTML code of the card.
    """
    ps = []
    if type(message) is list:
        for a_msg in message:
            ps.append(html.P(className='card-text', children=a_msg))
    else:
        ps = [html.P(className='card-text', children=message)]

    # return [html.Div(className='card', style={'min-width': '15rem', 'max-width': '25rem'}, children=[
    return [html.Div(className='card m-3', children=[
            html.Div(className='card-body', children=[
                html.H5(className='card-title', children=title)] + ps)])
        ]


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
    pol_dict = {1: 'single', 2: 'dual', 4: 'full'}
    rms = optimal_units(obs.thermal_noise(), [u.MJy, u.kJy, u.Jy, u.mJy, u.uJy])
    rms_time = optimal_units(obs.thermal_noise()/np.sqrt(obs.inttime/obs.duration),
                             [u.MJy, u.kJy, u.Jy, u.mJy, u.uJy])
    rms_channel = optimal_units(rms*np.sqrt(obs.subbands*obs.channels),
                                [u.MJy, u.kJy, u.Jy, u.mJy, u.uJy])

    ants_up = obs.is_visible()
    ant_no_obs = []
    for an_ant in ants_up:
        if len(ants_up[an_ant][0]) == 0:
            ant_no_obs.append(an_ant)

    ant_text = ', '.join(obs.stations.keys()) + '.'

    cards = []
    # The time card
    temp_msg = [f"{fx.print_obs_times(obs)}."]
    temp_msg += [f"Observing for {optimal_units(obs.duration, [u.h, u.min, u.s, u.ms]):.3n} in total, with {optimal_units(obs.ontarget_time, [u.h, u.min, u.s, u.ms]):.3n} on target."]
    temp_msg += [f"With a time integration of {optimal_units(obs.inttime, [u.s,u.ms,u.us]):.2n} the expected FITS file size is {optimal_units(obs.datasize(), [u.TB, u.GB, u.MB, u.kB]):.3n}."]
    cards += create_sensitivity_card('Observing Time', temp_msg)
    # Antennas
    longest_bl = obs.longest_baseline()[1]
    # Using dummy units to allow the conversion
    longest_bl_lambda = longest_bl/obs.wavelength
    longest_bl_lambda = optimal_units(longest_bl_lambda*u.m, [u.Gm, u.Mm, u.km])
    temp_msg = [f"{len(ants_up)-len(ant_no_obs)} participating antennas: {ant_text}"]
    if len(ant_no_obs) > 0:
        temp_msg += [html.P(className='text-danger', children=f"Note that {', '.join(ant_no_obs)} cannot observe the source.")]

    temp_msg += [f"The longest (projected) baseline is {optimal_units(longest_bl, [u.km, u.m]):.5n} ({longest_bl_lambda.value:.3n} {longest_bl_lambda.unit.name[0]}lambda)."]
    synthbeam = obs.synthesized_beam()
    synthbeam_units = optimal_units(synthbeam['bmaj'], [u.arcsec, u.mas, u.uas]).unit
    temp_msg += [html.P([f"The expected synthesized beam will be approx. {synthbeam['bmaj'].to(synthbeam_units).value:.2n} x {synthbeam['bmin'].to(synthbeam_units):.2n}", html.Sup("2"), \
            f", PA = {synthbeam['pa']:.3n}."])]
    cards += create_sensitivity_card('Antennas', temp_msg)

    # Frequency
    temp_msg = [f"Observing at a central frequency of {optimal_units(obs.frequency, [u.GHz, u.MHz]):.3n} ({optimal_units(obs.wavelength, [u.m, u.cm, u.mm]):.2n})."]
    temp_msg += [f"The total bandwidth of {optimal_units(obs.bandwidth, [u.GHz, u.MHz, u.kHz]):.3n} will be divided into {obs.subbands} subbands of {optimal_units(obs.bandwidth/obs.subbands, [u.GHz, u.MHz, u.kHz]):.3n} each, with {obs.channels} channels ({optimal_units(obs.bandwidth/(obs.subbands*obs.channels), [u.GHz, u.MHz, u.kHz, u.Hz]):.3n} wide)."]
    temp_msg += [f"Recording {pol_dict[obs.polarizations]} circular polarization."]
    cards += create_sensitivity_card('Frequency Setup', temp_msg)

    # RMS
    temp_msg = [f"Considering the sensitivities of the antennas, the estimated thermal noise is {rms:.3n}/beam."]
    temp_msg += [f"This would imply a rms of {rms_channel:.3n}/beam per spectral channel, or approx. {rms_time:.3n}/beam per time integration ({optimal_units(obs.inttime, [u.s,u.ms,u.us]):.3n})."]
    cards += create_sensitivity_card('Sensitivity', temp_msg)

    # resolution and FITS files
    # temp_msg = [f""]
    # cards += create_sensitivity_card('Resolution')

    # FoV smearing
    shortest_bl = obs.shortest_baseline()[1]
    largest_ang_scales = ((2.063e8*u.mas)*(obs.wavelength/shortest_bl)).to(u.mas)
    temp_msg = [f"The Field of View would be limited by time smearing to {optimal_units(obs.time_smearing(), [u.arcmin, u.arcsec]):.3n} and by frequency smearing to {optimal_units(obs.bandwidth_smearing(), [u.arcmin, u.arcsec]):.3n} (for a 10% loss)."]
    temp_msg += [f"Considering the shortest baseline in the array ({optimal_units(shortest_bl, [u.km, u.m]):.5n}), you will filter out emission on angular scales larger than {optimal_units(largest_ang_scales, [u.arcmin, u.arcsec, u.mas]):.3n}."]
    cards += create_sensitivity_card('FoV limitations', temp_msg)

    return [html.Div(className='card-deck col-12', children=cards)]



@app.callback(Output('onsourcetime-label', 'children'),
              [Input('onsourcetime', 'value')])
def update_onsourcetime_label(onsourcetime):
    return f"% of on-target time ({onsourcetime}%)"


@app.callback(Output('antennas-div', 'children'),
              [Input('band', 'value'), Input('array', 'value'), Input('e-EVN', 'value')])
def select_antennas(selected_band, selected_networks, is_eEVN):
    """Given a selected band and selected default networks, it selects the associated
    antennas from the antenna list.
    """
    selected_antennas = []
    if is_eEVN:
        selected_antennas = [ant for ant in default_arrays['e-EVN'] \
                             if (all_antennas[ant].has_band(selected_band) and \
                                (all_antennas[ant].network == 'EVN'))]

        return [html.Div([html.Br(),
                html.Label(html.H4(f"{sorted_networks[an_array]}")),
                html.Br(),
                dcc.Checklist(id=f"list_stations_{an_array}",
                    options=[{'label': s.name, 'value': s.codename,
                    'disabled': (not s.has_band(selected_band)) or \
                        (not s.codename in default_arrays['e-EVN'])}
                    for s in all_antennas if s.network == an_array],
                    value=selected_antennas if an_array=='EVN' else [],
                    className='antcheck', labelClassName='form-check-label',
                    inputClassName='form-check-input')]) for an_array in sorted_networks]
    else:
        for an_array in selected_networks:
            selected_antennas += [ant for ant in default_arrays[an_array] \
                                    if all_antennas[ant].has_band(selected_band)]

        return [html.Div([html.Br(),
                html.Label(html.H4(f"{sorted_networks[an_array]}")),
                html.Br(),
                dcc.Checklist(id=f"list_stations_{an_array}",
                    options=[{'label': s.name, 'value': s.codename,
                    'disabled': not s.has_band(selected_band)}
                    for s in all_antennas if s.network == an_array],
                    value=[s.codename for s in all_antennas \
                            if (s.codename in selected_antennas) and (s.network == an_array)],
                    className='antcheck', labelClassName='form-check-label',
                    inputClassName='form-check-input')]) for an_array in sorted_networks]


@app.callback([Output('error_starttime', 'children'),
               Output('error_endtime', 'children')],
              [Input('starttime', 'value'), Input('endtime', 'value')])
def check_obstime(starttime, endtime):
    if starttime != 'DD/MM/YYYY HH:MM':
        try:
            time0 = Time(datetime.datetime.strptime(starttime, '%d/%m/%Y %H:%M'),
                         format='datetime')
        except ValueError as e:
            return 'Incorrect format (dd/mm/YYYY HH:MM)', dash.no_update

    if endtime != 'DD/MM/YYYY HH:MM':
        try:
            time1 = Time(datetime.datetime.strptime(endtime, '%d/%m/%Y %H:%M'),
                         format='datetime')
        except ValueError as e:
            return dash.no_update, 'Incorrect format (dd/mm/YYYY HH:MM)'

    if ('time1' in locals()) and ('time0' in locals()):
        if (time1 - time0) > 5*u.d:
            return ["Please, put an observation shorter than 5 d"]*2
        elif (time0 - time1) >= 0*u.d:
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
               State('starttime', 'value'),
               State('endtime', 'value'),
               State('source', 'value'),
               State('onsourcetime', 'value'),
               State('datarate', 'value'),
               State('subbands', 'value'),
               State('channels', 'value'),
               State('pols', 'value'),
               State('inttime', 'value'),
               State('list_stations_EVN', 'value'),
               State('list_stations_eMERLIN', 'value'),
               State('list_stations_VLBA', 'value'),
               State('list_stations_LBA', 'value'),
               State('list_stations_KVN', 'value'),
               State('list_stations_Other', 'value')])
def compute_observation(n_clicks, band, starttime, endtime, source, onsourcetime, datarate,
                        subbands, channels, pols, inttime, *ants):
                       # subbands, channels, pols, inttime, ants_evn, ants_emerlin,
                       # ants_vlba, ants_lba, ants_kvn, ants_other):
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
        time0 = Time(datetime.datetime.strptime(starttime, '%d/%m/%Y %H:%M'),
                     format='datetime', scale='utc')
    except ValueError as e:
        return alert_message("Incorrect format for starttime."), \
               "First, set correctly an observation in the previous tab.", \
               dash.no_update, dash.no_update, dash.no_update, dash.no_update

    try:
        time1 = Time(datetime.datetime.strptime(endtime, '%d/%m/%Y %H:%M'),
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
    # obs_times = time0 + np.arange(0, (time1-time0).to(u.min).value, 15)*u.min
    all_selected_antennas = list(itertools.chain.from_iterable(ants))
    try:
        obs = observation.Observation(target=target_source, times=obs_times, band=band,
                      datarate=datarate, subbands=subbands, channels=channels,
                      polarizations=pols, inttime=inttime, ontarget=onsourcetime/100.0,
                      stations=get_selected_antennas(all_selected_antennas))
        sensitivity_results = update_sensitivity(obs)
    except observation.SourceNotVisible:
        return alert_message([
                    html.P(["Your source cannot be observed within the arranged observation.",
                    html.Br(),
                    "The source is not visible for any of the selected antennas " \
                    + "in the given observing time."]),
                    html.P("Modify the observing time of select a different array to observe" \
                    + " this source.")], title="Warning!"), \
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
                        'mode': 'lines', 'hovertemplate': "Elev: %{y:.2n}º",
                        'name': obs.stations[ant].name})

    data_fig.append({'x': obs.times.datetime, 'y': np.zeros_like(obs.times)+10,
                     'mode': 'lines', 'hoverinfo': 'skip', 'name': 'Elev. limit 10º',
                     'line': {'dash': 'dash', 'opacity': 0.5, 'color': 'gray'}})
    data_fig.append({'x': obs.gstimes.value, 'y': np.zeros_like(obs.times)+20, 'xaxis': 'x2',
                     'mode': 'lines', 'hoverinfo': 'skip', 'name': 'Elev. limit 20º',
                     'line': {'dash': 'dot', 'opacity': 0.5, 'color': 'gray'}})
    return {'data': data_fig,
            'layout': {'title': 'Source elevation during the observation',
                       'hovermode': 'closest',
                       'xaxis': {'title': 'Time (UTC)', 'showgrid': False,
                                 'ticks': 'inside', 'showline': True, 'mirror': False,
                                 'hovermode': 'closest', 'color': 'black'},
                       'xaxis2': {'title': {'text': 'Time (GST)', 'standoff': 0},
                                  'showgrid': False, 'overlaying': 'x',
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
                         # 'mode': 'markers', 'hoverinfo': "skip",
                         'mode': 'markers', 'marker': {'symbol': "41"}, 'hoverinfo': "skip",
                         'name': obs.stations[ant].name})

    return {'data': data_fig,
            'layout': {'title': 'Source visible (>10 deg) during the observation',
                       'xaxis': {'title': 'Time (UTC)', 'showgrid': False,
                                 'ticks': 'inside', 'showline': True, 'mirror': "all",
                                 'hovermode': 'closest', 'color': 'black'},
                       'yaxis': {'ticks': '', 'showline': True, 'mirror': True,
                                 'showticklabels': False, 'zeroline': False,
                                 'showgrid': False, 'hovermode': 'closest',
                                 'startline':False}}}




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


