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
import plotly.graph_objs as go
from datetime import datetime as dt
from astropy.time import Time
from astropy import coordinates as coord
from astropy import units as u
# Tweak to not let astroplan crashing...

from astropy.utils import iers

from astroplan import FixedTarget
from src import freqsetups as fs
from src import stations
from src import functions as fx
from src import observation


current_directory = path.dirname(path.realpath(__file__))
# stationList =  stations.Stations()
# stationList.add_from_file(current_directory+'/station_location.txt')

iers.conf.auto_download = False
iers.conf.auto_max_age = None
local_iers = f"{current_directory}/data/finals2000A.all"
iers.IERS.iers_table = iers.IERS_A.open(local_iers, cache=True)
iers.conf.remote_timeout = 300
# iers.IERS.iers_table = iers.IERS.open(cache=True)
# iers.IERS.iers_table = iers.IERS_A.open(iers.IERS_A_URL)
# iers.Conf.iers_auto_url.set('ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all')



all_antennas = fx.get_stations_from_file(f"{current_directory}/data/station_location.txt")
sorted_networks = ('EVN', 'eMERLIN', 'VLBA', 'LBA', 'KVN', 'Other')
default_arrays = {'EVN': ['Ef', 'Hh', 'Jb2', 'Mc', 'Nt', 'Ur', 'On', 'Sr', 'T6', 'Tr',
                          'Ys', 'Wb', 'Bd', 'Sv', 'Zc', 'Ir'],
          'e-EVN': ['Ef', 'Hh', 'Ir', 'Jb2', 'Mc', 'Nt', 'On', 'T6', 'Tr', 'Ys', 'Wb',
                    'Bd', 'Sv', 'Zc', 'Ir', 'Sr', 'Ur'],
          'eMERLIN': ['Cm', 'Kn', 'Pi', 'Da', 'De'],
          # 'LBA': ['ATCA', 'Pa', 'Mo', 'Ho', 'Cd', 'Td', 'Ww'],
          'LBA': ['ATCA', 'Pa', 'Mo', 'Ho', 'Cd', 'Td'],
          'VLBA': ['Br', 'Fd', 'Hh', 'Kp', 'La', 'Mk', 'Nl', 'Ov', 'Pt', 'Sc'],
          'KVN': ['Ky', 'Ku', 'Kt'],
          'Global VLBI': ['Ef', 'Hh', 'Jb2', 'Mc', 'Nt', 'Ur', 'On', 'Sr', 'T6',
                          'Tr', 'Ys', 'Wb', 'Bd', 'Sv', 'Zc', 'Ir', 'Br', 'Fd', 'Hn',
                          'Kp', 'La', 'Mk', 'Nl', 'Ov', 'Pt', 'Sc'],
          'GMVA': ['Ef', 'Mh', 'On', 'Ys', 'Pv', 'Br', 'Fd', 'Kp', 'La', 'Mk', 'Nl',
                   'Ov', 'Pt'],
          'EHT': ['ALMA', 'Pv', 'LMT', 'PdB', 'SMA', 'JCMT', 'APEX', 'SMTO', 'SPT']}

# Initial values
target_source = observation.Source('1h2m3s +50d40m30s', 'Source')
# obs_times = Time('1967-04-17 10:00') + np.arange(0, 600, 15)*u.min
selected_band = '18cm'

sensitivity_results_template = """

{band:.2n} ({freq:.2n}) observations with the following antennas:

{antennas}

{sb} ({sbbw:.2n}) subbands with {ch} channels each and {pols} polarization.

Total bandwidth of the observation: {bandwidth:.3n}
({bandwidth_channel:.3n} per spectral channel)


**Estimated image thermal noise** (assuming {ttarget:.3n} on target):  {noise:.3n}


Estimated rms thernal noise per spectral channel: {noise_channel:.3n}

**Resulting FITS file size**: {filesize:.3n}

(note that this is only an estimation)


**Smearing in the Field of View** (for a 10% loss):

Field of View limited by bandwidth-smearing to: {bw_smearing:.3n}

Field of View limited by time-smearing to: {t_smearing:.3n}

"""




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


#####################  This is the webpage layout
app.layout = html.Div([
    html.Div([
        html.H1('EVN Source Visibility'),
        # html.Img(src='http://www.ira.inaf.it/evnnews/archive/evn.gif')
    ]),
    # ], className='banner'),
    html.Div([html.Br()]), #style={'clear': 'both', 'margin-top': '20px'}),
    # First row containing all buttons/options, list of telescopes, and button with text output
    dcc.ConfirmDialog(id='global-error', message=''),
    # Elements in second column (checkboxes with all stations)
    html.Div(className='container-fluid', children=[
        dcc.Tabs([
            dcc.Tab(label='Observation Setup', className='container', children=[
                # Elements in first column ()
                html.Div(className='row-cols-3', children=[
                html.Div(className='col-sm-3', style={'max-width': '300px','float': 'left',
                                                      'padding': '2%'}, children=[
                    html.Div(className='form-group', children=[
                        html.Label('Observing Band'),
                        dcc.Dropdown(id='band', persistence=True,
                                 options=[{'label': fs.bands[b], 'value': b} for b \
                                in fs.bands], value='18cm'),
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('Select default VLBI Network(s)'),
                        dcc.Dropdown(id='array', options=[{'label': n, 'value': n} \
                                for n in default_arrays], value=['EVN'], multi=True),
                    ]),
                    html.Div(className='input-group-prepend', children=[
                        dcc.Checklist(id='e-EVN', className='checkbox', persistence=True,
                                      options=[{'label': ' e-EVN (real-time) mode?',
                                                'value': 'e-EVN'}], value=[]),
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
                        # dcc.Input(id='source', value='hh:mm:ss dd:mm:ss', type='text',
                        dcc.Input(id='source', value='01:20:00 +50:40:30', type='text',
                                  className='form-control', placeholder="hh:mm:ss dd:mm:ss",
                                  persistence=True),
                        html.Small(id='error_source', style={'color': 'red'},
                                   className='form-text text-muted'),
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label(id='onsourcetime-label',
                                   children='Percent. of on-target time'),
                        dcc.Slider(id='onsourcetime', min=20, max=100, step=5, value=75,
                                   marks= {i: str(i) for i in range(20, 101, 10)},
                                   persistence=True),
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('Datarate per station (in Mbps)'),
                        dcc.Dropdown(id='datarate', placeholder="Select a datarate...",
                                     options=[{'label': str(dr), 'value': dr} \
                                     for dr in fs.data_rates], value=1024, persistence=True),
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('Number of subbands'),
                        dcc.Dropdown(id='subbands', placeholder="Select no. subbands...",
                                     options=[{'label': str(sb), 'value': sb} \
                                     for sb in fs.subbands], value=8, persistence=True),
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('Number of spectral channels'),
                            dcc.Dropdown(id='channels', placeholder="Select no. channels...",
                                         options=[{'label': str(ch), 'value': ch} \
                                         for ch in fs.channels], value=32, persistence=True),
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('Number of polarizations'),
                        dcc.Dropdown(id='pols', placeholder="Select polarizations...",
                                     options=[{'label': fs.polarizations[p], 'value': p} \
                                     for p in fs.polarizations], value=4, persistence=True),
                    ]),
                    html.Div(className='form-group', children=[
                        html.Label('Integration time (s)'),
                        dcc.Dropdown(id='inttime', placeholder="Select integration time...",
                                     options=[{'label': fs.inttimes[it], 'value': it} \
                                     for it in fs.inttimes], value=2, persistence=True),
                    ])
                ]),
                # html.Div(style={'margin-top': '20px'}, children=[
                html.Div(className='col-lg-7', style={'float': 'left'}, children=[
                    html.Div(id='antennas-div', className='container', children=[
                    # List with all antennas
                        html.Div(className='antcheck', children=[html.Br(),
                            html.Label(html.H4(f"{an_array}")),
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
                html.Div(className='col-sm-2', style={'float': 'left'}, children=[
                    html.Button('Compute Observation', id='antenna-selection-button',
                                className='btn btn-primary btn-lg',
                                # style={'margin': '5px 5px 5px 5px'}),
                                style={'padding': '5px', 'margin-top': '50px'}),
                ])
                ])
            ]),

            dcc.Tab(label='Sensitivity', children=[
                html.Div(className='col-md-8', children=[
                # Sensitivity calculations
                    dcc.Markdown(id='sensitivity-output',
                                 children="Set the observation first.")
                ])
            ]),
            dcc.Tab(label='Plots', children=[
                html.Div(className='col-md-8', children=[
                    # Elevation VS time
                    html.Div([
                        dcc.Graph(id='fig-elev-time')
                    ]),
                    # Antenna VS time (who can observe)
                    html.Div([
                        dcc.Graph(id='fig-ant-time')
                    ])
                ])
            ]),
            dcc.Tab(label='Images', children=[
                #  Images
                html.Div(className='col-md-8', children=[
                    dcc.Markdown(children="""To be implemented.
                        The uv coverage and expected dirty images will go here.""")
                ])
            ]),
            dcc.Tab(label='Help', children=[
                #  Documentation
                html.Div(className='col-md-8', children=[
                    dcc.Markdown(children="""The Help/doc.
                        All explanations and technical details will go here.""")
                ])
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


def update_sensitivity(an_obs):
    """Given the observation, it sets the text about all properties of the observation.
    """
    rms = optimal_units(an_obs.thermal_noise(), [u.MJy, u.kJy, u.Jy, u.mJy, u.uJy])
    rms_channel = optimal_units(rms*np.sqrt(an_obs.subbands*an_obs.channels),
                                [u.MJy, u.kJy, u.Jy, u.mJy, u.uJy])

    ants_up = an_obs.is_visible()
    ant_no_obs = []
    for an_ant in ants_up:
        if len(ants_up[an_ant][0]) == 0:
            ant_no_obs.append(an_ant)

    antennas_text = ', '.join(an_obs.stations.keys())
    if len(ant_no_obs) > 0:
        antennas_text += f"\n (note that {', '.join(ant_no_obs)} cannot observe the source)."
    return sensitivity_results_template.format(
            band=optimal_units(an_obs.wavelength, [u.m, u.cm, u.mm]),
            freq=optimal_units(an_obs.frequency, [u.GHz, u.MHz]),
            antennas=antennas_text,
            sb=an_obs.subbands,
            sbbw=optimal_units(an_obs.bandwidth/an_obs.subbands, [u.GHz, u.MHz, u.kHz]),
            ch=an_obs.channels,
            pols={1: 'single', 2: 'dual', 4: 'full'}[an_obs.polarizations],
            ttarget=optimal_units(an_obs.ontarget_fraction*(an_obs.times[-1]-an_obs.times[0]),
                                  [u.h, u.min, u.s, u.ms]),
            noise=rms,
            noise_channel=rms_channel,
            filesize=optimal_units(an_obs.datasize(), [u.TB, u.GB, u.MB, u.kB]),
            bw_smearing=optimal_units(an_obs.bandwidth_smearing(), [u.arcmin, u.arcsec]),
            t_smearing=optimal_units(an_obs.time_smearing(), [u.arcmin, u.arcsec]),
            bandwidth=optimal_units(an_obs.bandwidth, [u.GHz, u.MHz, u.kHz]),
            bandwidth_channel=optimal_units(an_obs.bandwidth/(an_obs.subbands*an_obs.channels),
                                            [u.GHz, u.MHz, u.kHz, u.Hz]))



@app.callback(Output('onsourcetime-label', 'children'),
              [Input('onsourcetime', 'value')])
def update_onsourcetime_label(onsourcetime):
    return f"Percent. of on-target time ({onsourcetime}%)"


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
                html.Label(html.H4(f"{an_array}")),
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
                html.Label(html.H4(f"{an_array}")),
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
            return ["Please, put a time range smaller than 5 days."]*2

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



@app.callback([Output('sensitivity-output', 'children'),
               Output('fig-elev-time', 'figure'),
               Output('fig-ant-time', 'figure'), Output('global-error', 'message')],
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
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update
    try:
        target_source = observation.Source(convert_colon_coord(source), 'Source')
    except ValueError as e:
        return f"""Incorrect format for source coordinates:
        {source} found but 'hh:mm:ss dd:mm:ss' expected.
        """, dash.no_update, dash.no_update, dash.no_update
    try:
        time0 = Time(datetime.datetime.strptime(starttime, '%d/%m/%Y %H:%M'),
                     format='datetime')
    except ValueError as e:
        return "Incorrect format for starttime.", dash.no_update, dash.no_update, dash.no_update

    try:
        time1 = Time(datetime.datetime.strptime(endtime, '%d/%m/%Y %H:%M'),
                     format='datetime')
    except ValueError as e:
        return "Incorrect format for endtime.", dash.no_update, dash.no_update, dash.no_update

    if time0 >= time1:
        return "The start time of the observation must be earlier than the end time.", \
                dash.no_update, dash.no_update, dash.no_update

    if (time1 - time0) > 5*u.d:
        return "Please, put a time range smaller than 5 days.", \
                dash.no_update, dash.no_update, dash.no_update

    # try:
    # TODO: this should not be hardcoded...
    obs_times = time0 + np.linspace(0, (time1-time0).to(u.min).value, 50)*u.min
    # obs_times = time0 + np.arange(0, (time1-time0).to(u.min).value, 15)*u.min
    all_selected_antennas = list(itertools.chain.from_iterable(ants))
    obs = observation.Observation(target=target_source, times=obs_times, band=band,
                      datarate=datarate, subbands=subbands, channels=channels,
                      polarizations=pols, inttime=inttime, ontarget=onsourcetime/100.0,
                      stations=get_selected_antennas(all_selected_antennas))

    # except Exception as e:
    #     return dash.no_update, dash.no_update, dash.no_update, error_text(e)


    # return update_sensitivity(obs), dash.no_update, dash.no_update
    return update_sensitivity(obs), get_fig_ant_elev(obs), get_fig_ant_up(obs), dash.no_update



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
    data_fig.append({'x': obs.times.datetime, 'y': np.zeros_like(obs.times)+20,
                     'mode': 'lines', 'hoverinfo': 'skip', 'name': 'Elev. limit 20º',
                     'line': {'dash': 'dot', 'opacity': 0.5, 'color': 'gray'}})
    return {'data': data_fig,
            'layout': {'title': 'Source elevation during the observation',
                       'xaxis': {'title': 'Time (UTC)', 'showgrid': False,
                                 'ticks': 'inside', 'showline': True, 'mirror': "all",
                                 'hovermode': 'closest', 'color': 'black'},
                       'yaxis': {'title': 'Elevation (degrees)', 'range': [0., 92.],
                                 'ticks': 'inside', 'showline': True, 'mirror': "all",
                                 'showgrid': False, 'hovermode': 'closest'}}}



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
            'layout': {'title': 'Source visible during the observation',
                       'xaxis': {'title': 'Time (UTC)', 'showgrid': False,
                                 'ticks': 'inside', 'showline': True, 'mirror': "all",
                                 'hovermode': 'closest', 'color': 'black'},
                       'yaxis': {'ticks': '', 'showline': True, 'mirror': True,
                                 'showticklabels': False, 'zeroline': False,
                                 'showgrid': False, 'hovermode': 'closest',
                                 'startline':False}}}





# @app.callback(Output('my-graph', 'figure'), [Input('my-dropdown', 'value')])
# def update_graph(selected_dropdown_value):
#     #df = web.DataReader(
#     #    selected_dropdown_value, data_source='google',
#     #    start=dt(2017, 1, 1), end=dt.now())
#     return {
#         'data': [{
#             'x': [1, 2, 3, 4, 5, 6],
#             'y': [3, 4, 5, 4, 5, 6]
#         }]
#     }





if __name__ == '__main__':
    # app.run_server(host='0.0.0.0', debug=True)
    # app.run_server(debug=True)
    app.run_server(host='0.0.0.0')


