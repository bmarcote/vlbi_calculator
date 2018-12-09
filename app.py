#! /usr/bin/env python3

from os import path
from time import sleep
import numpy as np
import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
from datetime import datetime as dt
from astroplan import Observer
import stations
import util_functions as uf
import plot_functions as pf

current_directory = path.dirname(path.realpath(__file__))
stationList =  stations.Stations()
stationList.add_from_file(current_directory+'/station_location.txt')


arrays = {'EVN': ['EF', 'HH', 'JB2', 'MC', 'NT', 'UR', 'ON', 'SR', 'T6', 'TR', 'YS',
                  'WB', 'BD', 'SV', 'ZC', 'IR', 'MH'],
          'e-EVN': ['EF', 'HH', 'IR', 'JB2', 'MC', 'NT', 'ON', 'T6', 'TR', 'YS', 'WB'],
          'LBA': ['ATCA', 'PA', 'MO', 'HO', 'CD', 'TD70', 'WW'],
          'VLBA': ['VLBA-BR', 'VLBA-FD', 'VLBA-HN', 'VLBA-KP', 'VLBA-LA', 'VLBA-MK',
                   'VLBA-NL', 'VLBA-OV', 'VLBA-PT', 'VLBA-SC'],
          'Global': ['EF', 'HH', 'JB2', 'MC', 'NT', 'UR', 'ON', 'SR', 'T6', 'TR', 'YS',
                     'WB', 'BD', 'SV', 'ZC', 'IR', 'MH', 'VLBA-BR', 'VLBA-FD', 'VLBA-HN',
                     'VLBA-KP', 'VLBA-LA', 'VLBA-MK', 'VLBA-NL', 'VLBA-OV', 'VLBA-PT',
                     'VLBA-SC'],
          'GMVA': ['EF', 'MH', 'ON', 'YS', 'PV', 'VLBA-BR', 'VLBA-FD', 'VLBA-KP',
                   'VLBA-LA', 'VLBA-MK', 'VLBA-NL', 'VLBA-OV', 'VLBA-PT'],
          'EHT': ['ALMA', 'PV', 'LMT', 'PdB', 'SMA', 'JCMT', 'APEX', 'SMTO', 'SPT']}

bands = {'92cm': 'P - 92cm', '49cm': 'P - 49cm', '30cm': 'UHF - 30cm', '21cm': 'L - 21cm',
         '18cm': 'L - 18cm', '13cm': 'S - 13cm', '6cm': 'C - 6cm', '5cm': 'C - 5cm',
         '3.6cm': 'X - 3.6cm', '2cm': 'U - 2cm', '1.3cm': 'K - 1.3cm', '0.9cm': 'Ka - 0.9cm',
         '0.7cm': 'Q - 0.7cm', '0.3cm': 'W - 0.3cm', '0.1cm': '0.1cm'}

data_rates = [4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8]


external_stylesheets = ['http://jive.eu/~marcote/style.css']
n_timestamps = 70 # Number of points (timestamps) for the whole observations.





app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
# app = dash.Dash(__name__)

app.config.requests_pathname_prefix = ''
# app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})
app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/brPBPO.css"})


#####################  This is the webpage layout
app.layout = html.Div([
    html.Div([
        html.H1('EVN Calculator'),
        # html.Img(src='http://www.ira.inaf.it/evnnews/archive/evn.gif')
    ], className='banner'),
    html.Div([html.Br()], style={'clear': 'both', 'margin-top': '100px'}),
    # First row containing all buttons/options, list of telescopes, and button with text output
    html.Div([
        # Elements in first column ()
        html.Div([
            html.Label('Array'),
            html.Br(),
            dcc.Dropdown(id='array', options=[{'label': a, 'value': a} for a in arrays],
                         value='EVN'),# multi=True),
            html.Br(),
            html.Label('Observing Band'),
            html.Br(),
            dcc.Dropdown(id='band', options=[{'label': bands[b], 'value': b} for b in bands],
                         value='18cm'),
            html.Br(),
            html.Label('Data rate (in Mbps)'),
            html.Br(),
            dcc.Dropdown(id='data_rate', options=[{'label': str(dr), 'value': dr} for dr in data_rates],
                         value=1024),
            html.Br(),
            html.Label('Source Coordinates'),
            html.Br(),
            # dcc.Input(id='source', value='hh:mm:ss dd:mm:ss', type='text', style={'width': '80%'}),
            dcc.Input(id='source', value='00:00:00 00:00:00', type='text', className='input_text'),
            html.Br(),
            html.Label('Start of observation (UTC)'),
            html.Br(),
            # dcc.Input(id='starttime', value='DD/MM/YYYY HH:MM', type='text', style={'width': '80%'}),
            dcc.Input(id='starttime', value='01/04/2010 12:00', type='text', className='input_text'),
            html.Br(),
            html.Label('Duration (in hours)'),
            html.Br(),
            dcc.Input(id='duration', value='24.0', type='text', className='input_text'),
            html.Br()
        ], className='column left'),
        # Elements in second column (checkboxes with all stations)
        html.Div([ # IN OPTIONS DICT,  ADD 'DISABLED' : bool  FOR DISABLED CHECKBOX. checked=checked
            dcc.Checklist(id='stations_check', options=[{'label': s.name, 'value': s.code}
                                  for s in stationList.stations], values=[],
                                  className='check_stations')
        ], className='column middle'),
        html.Div([
            html.Button('Calculate', id='button'),
            dcc.Markdown(''' ''', id='output-button'),
        ], className='column right')
    ]),
    # Second row with Plots
    html.Div([
        html.Div([
            dcc.Graph(id='plot-elevation', className='plot'),
        ]),
        html.Div([
            dcc.Graph(id='plot-visibility', className='plot')
        ]),
    ], className='element')
])

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



@app.callback(
    dash.dependencies.Output('stations_check', 'values'),
    [dash.dependencies.Input('array', 'value'),
    dash.dependencies.Input('band', 'value')])
def select_antennas(array, band):
    """Once a button for a specific array is pressed, we select the corresponding
    antennas.
    """
    return [a for a in arrays[array] if stationList[a].has_frequency(band)]


@app.callback(
    dash.dependencies.Output('stations_check', 'options'),
    [dash.dependencies.Input('array', 'value'),
    dash.dependencies.Input('band', 'value')])
def disable_antennas(array, band):
    return [{'label': s.name, 'value': s.code, 'disabled': not s.has_frequency(band),
            'className': 'disabled' if not s.has_frequency(band) else 'enabled'} for s
            in stationList.stations]


@app.callback(
    dash.dependencies.Output('output-button', 'children'),
    [dash.dependencies.Input('button', 'n_clicks')],
    [dash.dependencies.State('stations_check', 'values'),
    dash.dependencies.State('band', 'value'),
    dash.dependencies.State('data_rate', 'value'),
    dash.dependencies.State('source', 'value'),
    dash.dependencies.State('starttime', 'value'),
    dash.dependencies.State('duration', 'value')])
def calculate_everything(button, array_codes, band, data_rate, source, starttime, duration):
    if duration == '0.0' or not (float(duration) > 0.0):
        return 'Set a duration and press "Calculate"'

    duration = float(duration)
    antennas = stationList.get_stations_with_codes(array_codes)
    # If source and/or start of observation are not set, then regular sensitivity calc.
    if source == 'hh:mm:ss dd:mm:ss' or starttime == 'DD/MM/YYYY HH:MM':
        obs_times = uf.get_obs_times(uf.get_time('01/01/1993 03:14'), duration, duration)
        ant_times = [obs_times]*len(antennas)
        # Dummy elevation (40 deg) for all stations
        elevations = np.ones((len(antennas), len(obs_times)))*40
    else:
        # NOTE: Check for source and startime format!!!!
        obs_times = uf.get_obs_times(uf.get_time(starttime), duration, duration/n_timestamps)
        ant_times, elevations = uf.times_elev_with_source_above(uf.get_coordinates(source), antennas, obs_times)

    # Remove all stations that do not observe the source
    antennas = [a for a,t in zip(antennas,ant_times) if len(t) > 0]
    if len(antennas) > 0:
        ant_times = [t for t in ant_times if len(t) > 0]
        elevations = [e for e in elevations if len(e) > 0]
        rms = uf.get_thermal_noise(antennas, ant_times, data_rate*1e6, band)
        return formatted_text(rms, uf.get_bandwidth_smearing(), uf.get_time_smearing())


def formatted_text(rms, fov_bandwidth=None, fov_time=None):
    # rms is passed as Jy/beam. Here I pick the best units
    if rms > 0.99:
        rms_units = 'Jy/beam'
    elif rms > 0.00099:
        rms_units = 'mJy/beam'
        rms *= 1000
    else:
        rms_units = '&mu;Jy/beam'
        rms *= 1e6

    noise_str = '''
    ### Thermal noise

    The image thermal noise is estimated to be {0:.2f} {1} (at 1-sigma level) using natural weighting.

    '''
    band_str = '''
    Due to bandwidth smearing the FoV is limited to {} arcsec.

    '''
    time_str = '''
    Due to time smearing the FoV is limited to {} arcsec.

    '''
    fov = '''
    ### Limitations in field of view (FoV)

    {0}
    {1}

    These values have been calculated for 10% loss in the response of a point source, and they give
    the FoV radius from the pointing center.
    '''

    if (fov_bandwidth is None) and (fov_time is None):
        return noise_str.format(rms, rms_units).replace('  ', '')
    else:
        fov_f = fov.format(band_str.format('{:.2f}'.format(fov_bandwidth) if fov_bandwidth is not None else ''),
                           time_str.format('{:.2f}'.format(fov_time) if fov_time is not None else ''))
        return (noise_str.format(rms, rms_units) + fov_f).replace('  ', '')


@app.callback(
    dash.dependencies.Output('plot-elevation', 'figure'),
    [dash.dependencies.Input('button', 'n_clicks')],
    [dash.dependencies.State('stations_check', 'values'),
    dash.dependencies.State('band', 'value'),
    dash.dependencies.State('source', 'value'),
    dash.dependencies.State('starttime', 'value'),
    dash.dependencies.State('duration', 'value')])
def update_plot_elevation(button, antenna_codes, bands, source, starttime, duration):
    if (source == 'hh:mm:ss dd:mm:ss' or starttime == 'DD/MM/YYYY HH:MM' or float(duration) <= 0.0):
        return None

    obs_times = uf.get_obs_times(uf.get_time(starttime), float(duration),
                                 float(duration)/n_timestamps)
    antennas = stationList.get_stations_with_codes(antenna_codes)
    ant_times, elevations = uf.times_elev_with_source_above(uf.get_coordinates(source), antennas, obs_times)
    antennas = [a for a,t in zip(antennas,ant_times) if len(t) > 0]
    if len(antennas) > 0:
        ant_times = [t for t in ant_times if len(t) > 0]
        elevations = [e for e in elevations if len(e) > 0]
        return pf.time_elevation(ant_times, elevations, [s.name for s in antennas], xlim=(obs_times[0].datetime, obs_times[-1].datetime))


@app.callback(
    dash.dependencies.Output('plot-visibility', 'figure'),
    [dash.dependencies.Input('button', 'n_clicks')],
    [dash.dependencies.State('stations_check', 'values'),
    dash.dependencies.State('band', 'value'),
    dash.dependencies.State('source', 'value'),
    dash.dependencies.State('starttime', 'value'),
    dash.dependencies.State('duration', 'value')])
def update_plot_visibility(button, antenna_codes, bands, source, starttime, duration):
    if (source == 'hh:mm:ss dd:mm:ss' or starttime == 'DD/MM/YYYY HH:MM' or float(duration) <= 0.0):
        return None

    obs_times = uf.get_obs_times(uf.get_time(starttime), float(duration),
                                 float(duration)/n_timestamps)
    antennas = stationList.get_stations_with_codes(antenna_codes)
    ant_times, elevations = uf.times_elev_with_source_above(uf.get_coordinates(source), antennas, obs_times)
    antennas = [a for a,t in zip(antennas,ant_times) if len(t) > 0]
    if len(antennas) > 0:
        ant_times = [t for t in ant_times if len(t) > 0]
        elevations = [e for e in elevations if len(e) > 0]
        return pf.time_visibility(ant_times, elevations, [s.name for s in antennas], xlim=(obs_times[0].datetime, obs_times[-1].datetime))



if __name__ == '__main__':
    app.run_server(host='0.0.0.0', debug=True)
    # app.run_server(host='0.0.0.0')


