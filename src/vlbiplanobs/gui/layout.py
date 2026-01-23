from dash import html, dcc
from dash_bootstrap_components import Modal
from vlbiplanobs import freqsetups as fs
from vlbiplanobs.gui import inputs, outputs


def top_banner(app) -> html.Div:
    # return html.Div(id='banner',
    return inputs.card([html.Div(className='card-body m-0 p-0', children=[
                       html.A(className='d-inline-block mr-md-auto',
                              href="https://www.evlbi.org",
                              children=[html.Img(height='70px',
                                                 src=app.get_asset_url("logo_evn.png"),
                                                 alt='European VLBI Network (EVN)',
                                                 className="d-inline-block align-middle")]),
                       html.Div(html.H2('EVN Observation Planner'),
                                        #html.Span(' (v4.6 beta)', style={'color': '#a01d26'})],
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
            inputs.card(inputs.duration()),
            inputs.card(inputs.networks(app)),
            inputs.card(inputs.antenna_list(app)),
            inputs.card(inputs.source_and_epoch_selection()),
            inputs.card(inputs.correlations())])


def compute_buttons(app) -> html.Div:
    return html.Div([html.Div(className='m-0 p-0', children=[
        html.Div(className='row d-flex m-0 p-0', children=[
            inputs.compute_button(),
            outputs.download_button_div(),
            dcc.Download(id="download-data")])])])


def outputs_column(app) -> html.Div:
    return html.Div(children=[
        # First just the warnings
        html.Div(id='user-message', className='m-0 p-0',
                 children=html.Blockquote(className='text-secondary text-bold ms-2 px-2',
                                          children=["Set your VLBI observation on the left and the press ",
                                                    html.Em('CALCULATE!'),
                                                    html.Footer(className='text-sm pt-2',
                                                                children="The details of the expected outcome will "
                                                                         "appear here")])),

        html.A(id='download-link'),
        html.Div(id='out-sun', className='m-0 p-0'),
        html.Div(id='out-phaseref', className='m-0 p-0'),
        html.Div(id='out-ant', className='m-0 p-0'),
        # Now the actual results
        html.Div(className='col-12 m-0 p-0', children=[
            html.Div(className='row d-flex m-0 p-0', children=[
                html.Div(outputs.obs_time(), className='col-6 m-0 p-0 d-flex align-items-stretch', id='div-card-time'),
                html.Div(outputs.summary_freq_res(), className='col-6 m-0 p-0 d-flex align-items-stretch', id='div-card-vel')]),
            html.Div(className='row d-flex m-0 p-0', children=[
                html.Div(Modal(id="sensitivity-baseline-modal", size='xl', is_open=False)),
                html.Div(outputs.rms(), className='col-6 m-0 p-0 d-flex align-items-stretch', id='card-rms'),
                html.Div(outputs.resolution(), className='col-6 m-0 p-0 d-flex align-items-stretch', id='card-resolution')])]),
        # style={'align-items': 'stretch'}),
        html.Div(outputs.plot_elevations(), id='out-elevations', className='m-0 p-0', hidden=True),
        html.Div(outputs.plot_uv_coverage(), id='out-uv-coverage', className='m-0 p-0', hidden=True),
        html.Div(className='col-12 m-0 p-0', children=html.Div(className='row d-flex m-0 p-0', children=[
            html.Div(outputs.field_of_view(), className='col-6 m-0 p-0 d-flex align-items-stretch', id='div-card-fov'),
            html.Div(outputs.data_size(), id='card-datasize', className='col-6 m-0 p-0 d-flex align-items-stretch')])),
        html.Div(outputs.plot_worldmap(), id='out-worldmap', className='m-0 p-0', hidden=True)
        # html.Div(className='col-12 m-0 p-0', children=html.Div(className='row d-flex m-0 p-0', children=[
        #     ]))
    ])
