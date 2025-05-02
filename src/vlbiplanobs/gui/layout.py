import functools
from typing import Optional, Union, Sequence
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
from vlbiplanobs.gui import inputs, outputs


def top_banner(app) -> html.Div:
    return html.Div(id='banner',
                    className='d-flex p-0 shadow-sm bg-white card m-2',
                    children=[html.Div(className='card-body m-0 p-3', children=[
                       html.A(className='d-inline-block mr-md-auto',
                              href="https://www.evlbi.org",
                              children=[html.Img(height='70px',
                                                 src=app.get_asset_url("logo_evn.png"),
                                                 alt='European VLBI Network (EVN)',
                                                 className="d-inline-block align-middle")]),
                       html.Div(html.H2('EVN Observation Planner'),
                                className='d-inline-block align-middle mx-auto mb-0 font-weight-bolder',
                                style={'text-align': 'center', 'margin': '0', "flex": "1"}),
                       html.A(className='d-inline-block ml-auto pull-right',
                              href="https://www.jive.eu",
                              children=[html.Img(src=app.get_asset_url("logo_jive.png"),
                                                 height='70px',
                                                 alt='Joint Institute for VLBI ERIC (JIVE)')])],
                                       style={'display': 'flex', 'align-items': 'center',
                                              'justify-content': 'space-between'})])


def inputs_column(app) -> html.Div:
    return html.Div(children=[
            inputs.card(inputs.pick_band(fs.bands)),
            inputs.card(inputs.networks(app)),
            inputs.card(inputs.antenna_list(app)),
            html.Div(className='col-12', children=html.Div(className='row d-flex m-0 p-2', children=[
                 inputs.card(inputs.source_selection(), classNameDiv='col-6 m-0 p-0'),
                 inputs.card(inputs.epoch_selection(), classNameDiv='col-6 m-0 p-0')])),
            inputs.card(inputs.correlations())])


def compute_buttons(app) -> html.Div:
    return html.Div(className='col-12', children=[
        html.Div(className='row d-flex', children=[
            inputs.compute_button(),
            outputs.download_button(),
            dcc.Download(id="download-data")])])


def outputs_column(app) -> html.Div:
    return html.Div(children=[
        # First just the warnings
        html.Div(id='user-message', className='col-12',
                 children=outputs.info_card(["Set your VLBI observation and press ", html.Em('CALCULATE!')],
                                            "The details of the expected outcome will appear here.")),
        html.Div(id='out-sun', className='col-12'),
        html.Div(id='out-phaseref', className='col-12'),
        html.Div(id='out-ant', className='col-12'),
        # Now the actual results
        html.Div(className='col-12 m-0 p-0', children=html.Div(className='row d-flex m-0 p-2', children=[
            html.Div(dbc.Modal(id="sensitivity-baseline-modal", size='xl', is_open=False)),
            html.Div(outputs.rms(), className='col-6 m-0 p-0', id='card-rms'),
            html.Div(outputs.resolution(), className='col-6 m-0 p-0', id='card-resolution')])),
        # style={'align-items': 'stretch'}),
        html.Div(outputs.plot_elevations(), id='out-elevations', className='col-12', hidden=True),
        html.Div(outputs.plot_uv_coverage(), id='out-uv-coverage', className='col-12', hidden=True),
        html.Div(className='col-12 m-0 p-0', children=html.Div(className='row d-flex m-0 p-2', children=[
                html.Div(outputs.field_of_view(), className='col-6 m-0 p-0', id='div-card-fov'),
                html.Div(outputs.summary_freq_res(), className='col-6 m-0 p-0', id='div-card-vel')])),
        html.Div(outputs.plot_worldmap(), id='out-worldmap', className='col-12', hidden=True)
    ])
