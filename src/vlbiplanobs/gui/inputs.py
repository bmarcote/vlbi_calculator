# import numpy as np
# from datetime import datetime as dt
# import enum
from typing import Optional, Union
from datetime import datetime as dt
# from astropy import units as u
# from fpdf import FPDF
from dash import Dash, html, dcc, callback, Output, Input, State
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
# import plotly.express as px
from vlbiplanobs import freqsetups as fs
from vlbiplanobs import stations
from vlbiplanobs import observation





def top_banner(app):
    return [html.Div(id='banner',
                     className='d-flex p-0 shadow-sm bg-white card',
                     children=[html.Div(className='card-body m-0 p-3', children=[
                        html.A(className='d-inline-block mr-md-auto',
                               href="https://www.evlbi.org",
                               children=[html.Img(height='70px',
                                                  src=app.get_asset_url("logo_evn.png"),
                                                  alt='European VLBI Network (EVN)',
                                                  className="d-inline-block align-middle"),
                                ]),
                        html.Div(html.H2('EVN Observation Planner'),
                                className='d-inline-block align-middle mx-auto mb-0 font-weight-bolder',
                                 style={'text-align': 'center', 'margin': '0', "flex": "1"}),
                        html.A(className='d-inline-block ml-auto pull-right',
                               href="https://www.jive.eu",
                               children=[html.Img(src=app.get_asset_url("logo_jive.png"),
                                                  height='70px',
                                                  alt='Joint Institute for VLBI ERIC (JIVE)')
                               ])
                     ], style={'display': 'flex', 'align-items': 'center',
                               'justify-content': 'space-between'})
                     ])
            ]


def card(children: Optional[list] = None, className: str = '') -> html.Div:
    return html.Div(className=' '.join(['card m-2', className]), #style={'margin': '10px'},
                    children=[html.Div(className='card-body', children=children)])


def compute_button() -> list:
    """Returns the button to compute the observation
    """
    return html.Div(html.Button('CALCULATE', id='compute-observation',
                                className='btn bg-gradient-info btn-lg mx-auto w-75 m-4 p-2 active',
                                style={'position': 'sticky', 'top': '20px', 'background-color': '#9DB7C4'}),
                    className='text-center', style={'position': 'relative'})

def results() -> list:
    return []

def pick_band(bands: dict[str, str]) -> list:
    """Returns the card allowing the user to pick up the band.
    """
    top_labels = [{'label': 'ν (GHz)', 'style': {'font-weight': 'bold'}}] + \
                 [b.split('or')[1].replace('GHz', '').strip() for b in bands.values()]
    bottom_labels = [{'label': 'λ (cm)', 'style': {'font-weight': 'bold'}}] + \
                    [b.split('or')[0].replace('cm', '').strip() for b in bands.values()]
    return [html.Div(className='row col-12',
                     children=[html.Div(className='col-12', children=[
                        html.Div(className='row align-items-bottom', children=[
                                 html.Div(className='col-8', children=[
                                   html.H4("Observing band", className='text-dark font-weight-bold mb-0'),
                                 ]),
                                 html.Div(className='col-4 text-right', children=[
                                    dbc.Switch(label='Show wavelenths', value=False, id='switch-band-label',
                                               persistence=True),
                                 ])
                        ]), html.Br(),
                        html.Div(dcc.Slider(min=0, max=len(top_labels), step=1, value=0,
                                  marks={i: l for i, l in enumerate(bottom_labels)},
                                  included=False, id='band-slider', persistence=True)),
                    ])])]


def network_entry(app, network: str, bands: dict[str, str]) -> html.Div:
    """Returns the DIV that shows a given VLBI network to be selected
    """
    if len(observation._NETWORKS[network].observing_bands) > 1:
        if len(observation._NETWORKS[network].observing_bands) > 2:
            temp = ', and'
        else:
            temp = ' and'
        wavelengths = ', '.join([b.replace('cm', '').strip() \
                                 for b in observation._NETWORKS[network].observing_bands[:-1]]) + \
                      f"{temp} {observation._NETWORKS[network].observing_bands[-1].replace('cm',
                                                                                           ' cm.').strip()}"
        frequencies = ', '.join([bands[b].split('or')[1].replace('GHz', '').strip() \
                                 for b in observation._NETWORKS[network].observing_bands[:-1]]) + temp + \
                      bands[observation._NETWORKS[network].observing_bands[-1]].split('or')[1] + '.'
    elif len(observation._NETWORKS[network].observing_bands) == 1:
        wavelengths = observation._NETWORKS[network].observing_bands[0]
        frequencies = bands[observation._NETWORKS[network].observing_bands[0]].split('or')[1].replace('GHz', '')
    else:
        wavelengths = 'N/A'
        frequencies = 'N/A'

    has_full_name = stations.Stations.get_network_full_name(network) != network
    return html.Div([dbc.Card([dbc.CardImg(src=app.get_asset_url(f"network-{network.lower()}.png"),
                                           top=False,
                                           style={'opacity': 1.0, 'object-fit': 'cover', 'height': '10rem'},
                                           className='full-background'),
                    dbc.CardImgOverlay(dbc.CardBody([html.Div(className='d-flex align-items-center',
                                                              children=[
                                                      html.H4(network,
                                                              className='text-bold text-white mb-0 '
                                                                        'mt-0 flex-grow-1'),
                                                      dbc.Switch(label='', value=False,
                                                                 id=f"network-{network}",
                                                                 className='form-check ml-auto',
                                                                 persistence=True)]),
                                                     #            className='form-check'),
                                                     # html.H3(network, className='text-bold mb-0'),
                                                     # TODO: adapt text size to width of the card
                                                     html.H6(stations.Stations.get_network_full_name(network),
                                                             className='mt-0 mb-0 text-white',
                                                             style={'white-space': 'nowrap',
                                                                    'text-overflow': 'ellipsis',
                                                                    'overflow': 'hidden'}) \
                                                             if has_full_name else '',
                                                     html.Label(f"Observes at {wavelengths}", hidden=True,
                                                                className='card-text text-white',
                                                                id=f"network-{network}-label-wav",
                                                                style={'line-height': '1.2'}),
                                                     html.Label(f"Observes at {frequencies}", hidden=False,
                                                                className='card-text text-white',
                                                                id=f"network-{network}-label-freq",
                                                                style={'line-height': '1.2'})],
                                                    # className='card-body text-start p-0 pt-0 w-100'),),
                                                    className='card-body text-start p-0 pt-0 w-100'),),
                     ], className='m-2 card-background', style={'min-width': '11rem', 'height': '8rem',
                                                                'overflow': 'hidden',
                                                'opacity': 1.0},
                     id=f"network-{network}-card"),],
                    className='col-4')
    # return html.Div(className='card card-background card-background-mask-secondary', id=f"network-{network}-div",
    #          children=[html.Div(#className='full-background',
    #                             style={'background-image': f"url('./assets/network-{network.lower()}.png')"},
    #                             children=[
    #                                 html.Div(className='card-body text-start p-3 w-100',
    #                                          children=[html.Div(className='mb-3 d-flex col-1',
    #                                                             children=[dbc.Switch(label='', value=False,
    #                                                                                  id=f"network-{network}")]
    #                                                             ),
    #                                                    html.Div(className='mb-3 d-flex col-11',
    #                                                             children=[html.H5(network,
    #                                                                               className='text-light')])])])



def networks(app, networks, bands) -> list:
    return [card([html.H4("Select default VLBI Network(s)", className='text-dark font-weight-bold mb-0 pl-2 ml-4'),
            html.Div([network_entry(app, network, bands) for network in networks], className='d-flex flex-wrap')])]


def antenna_list(app, antennas) -> list:
    """Returns the DIV that shows all the antennas that can be selected, and allows searching through them
    """
    # TODO: change collapdsed to true
    return [card(dbc.Accordion(start_collapsed=False, className='accordion-contents-open-ids', flush=True,
                               persistence=True, id='accordion-state', children=[
            dbc.AccordionItem(title=html.H4(['Manual Selection of Antennas   ',
                                             html.I(className='fa fa-solid fa-angle-down', id='accordion-ant')],
                                            className='text-dark font-weight-bold mb-0 accordion-header'),
                              children=[
                html.Div(id='antennas-div', className='container mb-3',
                         style={'display': 'inline-block', 'gap': '10px', 'flex-wrap': 'wrap'},
                         children=[#html.Div(style={'display': 'flex', 'align-items': 'center', 'width': '10rem', 'box-sizing': 'border-box', 'width': 'var(--switch-width)'}, children=
                                dbc.Checklist(value=[], id='switches-antennas', inline=True,
                                              persistence=True, switch=True,
                                              style={'display': 'grid',
                                                     'grid-template-columns': 'repeat(auto-fit, minmax(10rem, 1fr))'},
                                              options=[{'label': s.name, 'value': s.codename} for s in antennas])
            ])
        ])]))]


def source_selection() -> html.Div:
    return card([html.H4("Source To Observe", className='text-dark font-weight-bold mb-0'),
        dbc.Switch(label='Specify name or coordinates', value=False,
                   id='switch-specify-source', persistence=True),
        html.Div(id='source-selection-div', hidden=True, children=[
            html.Div(className='row', children=[
                dcc.Input(id='source-input', value=None, type='text',
                          className='form-control', placeholder="hh:mm:ss dd:mm:ss", persistence=True),
                html.Small(id='error_source', className='form-text text-muted')
            ]),
            html.Div(className='form-group', children=[
                html.Label('Percentage of observing time spent on target'),
                dcc.Slider(id='onsourcetime', min=20, max=100, step=10, value=70,
                           marks= {i: str(i) for i in range(20, 101, 10)},
                           persistence=True)
            ])
        ])
    ])


def epoch_selection() -> html.Div:
    return card([html.H4("Define Observing Epoch", className='text-dark font-weight-bold mb-0'),
        dbc.Switch(label='Specify epoch', value=False,
                   id='switch-specify-epoch', persistence=True),
        html.Div(id='epoch-selection-div', hidden=True, children=[
            html.Div(className='row', children=[
                html.Label('Start of observation (UTC)'),
                html.Div(className='row mx-0', children=[
                    html.Div(className='col-12 px-0 mx-0', children=[
                    # dcc.DatePickerSingle(id='starttime', date=None, min_date_allowed=dt(1900, 1, 1),
                    #                      max_date_allowed=dt(2100, 1, 1),
                    #                      display_format='DD-MM-YYYY (DDD)',
                    #                      placeholder='Start date',
                    #                      first_day_of_week=1,
                    #                      initial_visible_month=dt.today(),
                    #                      persistence=True,
                    #                      className='form-picker',
                    #                      style={'width': '100%'}),
                    dmc.MantineProvider(
                        dmc.DateTimePicker(id='starttime', className='form-picker',
                                           value=None, minDate=dt(1900, 1, 1),
                                           maxDate=dt(2100, 1, 1),
                                           valueFormat='DD-MM-YYYY hh:mm',
                                           placeholder='Start time', withSeconds=False,
                                           persistence=True, clearable=True,
                                           style={'width': '100%'})),
                    # dmc.DateInput(id='starttime', value=None, minDate=dt(1900, 1, 1),
                    #                      maxDate=dt(2100, 1, 1),
                    #                      valueFormat='DD-MM-YYYY',
                    #                      placeholder='Start date',
                    #                      # first_day_of_week=1,
                    #                      # initial_visible_month=dt.today(),
                    #                      persistence=True, clearable=True,
                    #                      className='form-picker',
                    #                      style={'width': '100%'})),
                #     ]), html.Div(className='col-6 mx-0 px-0', children=[
                #         dcc.Dropdown(id='starthour', placeholder="Start time (UTC)", value=None,
                #                      options=[{'label': f"{hm//60:02n}:{hm % 60:02n}", \
                #                                'value': f"{hm//60:02n}:{hm % 60:02n}"} \
                #                               for hm in range(0, 24*60, 15)],
                #                      persistence=True, className='form-hour col-12',
                #                      style={'width': '100%'},),
                    ]),
                ]),
                html.Div(className='row', children=[
                    html.Small(id='error_starttime', style={'color': 'red'},
                               className='form-text text-muted'),
                ])
            ])
        ]),
        html.Div(className='form-group', children=[
            dcc.Input(id='duration', value=None, type='number', className='form-control',
                      placeholder="Duration of the observation (in hours)",
                      persistence=True, inputMode='numeric', min=0.5, max=50),
            html.Small(id='error_duration', className='form-text text-danger')
        ])
    ])


def correlations() -> html.Div:
    """Creates the inputs to specify if the user wants continuum correlation and/or spectral line.
    """
    return card([html.H4("Observation Setup", className='text-dark font-weight-bold mb-0'),
                 html.Div(className='col-12', children=[html.Div(className='row d-flex align-items-bottom',
                                                                 children=[
                    html.Div(className='col-6', children=[
                        dbc.Switch(label='Real-time (e-EVN observation)', value=False,
                                   id='switch-specify-e-evn', persistence=True),
                        dbc.Switch(label='Optimize for spectral line', value=False,
                                   id='switch-specify-continuum', persistence=True),
                        dcc.Dropdown(id='datarate', placeholder="Maximum observing data rate...",
                                     options=[{'label': fs.data_rates[dr], 'value': dr} \
                                     for dr in fs.data_rates][::-1], value=2048, persistence=True),
                        html.Label(id='bandwidth-label', style={'color': '#999999'}, children='',
                                   htmlFor="datarate")
                    ]),
                    html.Div(className='col-6', children=[
                        html.Div(className='form-group', children=[
                            dcc.Dropdown(id='subbands', placeholder="Select no. subbands...",
                                         #className='form-control',
                                         options=[{'label': fs.subbands[sb], 'value': sb} \
                                         for sb in fs.subbands][::-1], value=8, persistence=True),
                            dcc.Dropdown(id='channels', placeholder="Select no. channels per subband...",
                                         options=[{'label': fs.channels[ch],
                                                   'value': ch} \
                                         for ch in fs.channels][::-1], value=64, persistence=True),
                            dcc.Dropdown(id='pols', placeholder="Polarizations...",
                                         options=[{'label': fs.polarizations[p], 'value': p} \
                                         for p in fs.polarizations], value=4, persistence=True),
                            dcc.Dropdown(id='inttime', placeholder="Integration time...",
                                         options=[{'label': fs.inttimes[it], 'value': it} \
                                         for it in fs.inttimes], value=2, persistence=True),
                        ]),
                    ])
                    ]),
                 ])])

        # html.Div(className='form-group', children=[
        #             dcc.Dropdown(id='datarate', placeholder="Maximum observing data rate...",
        #                          options=[{'label': fs.data_rates[dr], 'value': dr} \
        #                          for dr in fs.data_rates][::-1], value=2048, persistence=True),
        #             html.Label(id='bandwidth-label', style={'color': '#999999'}, children='',
        #                        htmlFor="datarate"),
        #             dcc.Dropdown(id='subbands', placeholder="Select no. subbands...", #className='form-control',
        #                          options=[{'label': fs.subbands[sb], 'value': sb} \
        #                          for sb in fs.subbands][::-1], value=8, persistence=True),
        #             dcc.Dropdown(id='channels', placeholder="Select no. channels per subband...",
        #                          options=[{'label': fs.channels[ch],
        #                                    'value': ch} \
        #                          for ch in fs.channels][::-1], value=64, persistence=True),
        #             dcc.Dropdown(id='pols', placeholder="Polarizations...",
        #                          options=[{'label': fs.polarizations[p], 'value': p} \
        #                          for p in fs.polarizations], value=4, persistence=True),
        #             dcc.Dropdown(id='inttime', placeholder="Integration time...",
        #                          options=[{'label': fs.inttimes[it], 'value': it} \
        #                          for it in fs.inttimes], value=2, persistence=True),
        #         ]),
        # ])
