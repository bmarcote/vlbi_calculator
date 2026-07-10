from dash import html, dcc
import dash_bootstrap_components as dbc
from vlbiplanobs import freqsetups as fs
from vlbiplanobs.gui import inputs, outputs


def top_banner(app) -> html.Div:
    """Return the top banner with EVN and JIVE logos.

    Parameters
    ----------
    app : Dash app
        Dash application instance.

    Returns
    -------
    html.Div
        Top banner component.
    """
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
                                                 alt='Joint Institute for VLBI ERIC (JIVE)',
                                                 className='logo-light'),
                                        html.Img(src=app.get_asset_url("logo_jive_white.png"),
                                                 height='70px',
                                                 alt='Joint Institute for VLBI ERIC (JIVE)',
                                                 className='logo-dark')])],
                                       style={'display': 'flex', 'align-items': 'center',
                                              'justify-content': 'space-between'})])


def inputs_column(app) -> html.Div:
    """Return the inputs column with all input cards.

    Parameters
    ----------
    app : Dash app
        Dash application instance.

    Returns
    -------
    html.Div
        Inputs column component.
    """
    return html.Div(children=[
            inputs.card(inputs.pick_band(fs.bands)),
            inputs.card(inputs.duration()),
            inputs.card(inputs.source_and_epoch_selection()),
            inputs.card(inputs.networks(app)),
            inputs.card(inputs.antenna_list(app)),
            inputs.card(inputs.correlations())])


def export_button_div() -> html.Div:
    """Return the Dash elements used by the export/import to/from Polaris.
    The Location is used to parse import parameters.
    The Alert pops up when the Button is clicked to show whether the export
    was done directly to Polaris or to clipboard.

    Returns
    -------
    html.Div
        Export/import components.
    """
    return html.Div(className='m-0 p-0', children=[
        dcc.Location(id='url', refresh=False),
        inputs.export_button(),
        dbc.Alert("Nothing to see here", id='export-alert',
              is_open=False, color='success',
              duration=5000, dismissable=True)
    ])

def compute_buttons(app) -> html.Div:
    """Return the compute buttons section.

    Parameters
    ----------
    app : Dash app
        Dash application instance.

    Returns
    -------
    html.Div
        Compute buttons component.
    """
    return html.Div([html.Div(className='m-0 p-0', children=[
        html.Div(className='row d-flex m-0 p-0', children=[
            inputs.compute_button(),
            outputs.download_button_div(),
            dcc.Download(id="download-data")])
        ])
    ])

def compute_buttons_realtime(app) -> html.Div:
    """Return the layout for real-time mode with loading spinner.

    The PDF download button is rendered per-target inside each output tab.

    Parameters
    ----------
    app : Dash app
        Dash application instance.

    Returns
    -------
    html.Div
        Loading spinner component.
    """
    return html.Div(className='m-0 p-0', children=[
        html.Div(className='row d-flex m-0 p-0 justify-content-center', children=[
            html.Div(className='col-12 text-center', children=[
                dbc.Spinner(id='loading', color='#a01d26',
                            children=html.Div(id='loading-div'))])])])


def outputs_column(app) -> html.Div:
    """Return the outputs column with user message and results container.

    Parameters
    ----------
    app : Dash app
        Dash application instance.

    Returns
    -------
    html.Div
        Outputs column component.
    """
    return html.Div(children=[
        # Top-level user message (shared across tabs)
        html.Div(id='user-message', className='m-0 p-0',
                 children=html.Blockquote(className='text-secondary text-bold ms-2 px-2',
                                          children=["Set your VLBI observation on the left to see the "
                                                    "details of the expected outcome here.",
                                                    html.Footer(className='text-sm pt-2',
                                                                children="Add target sources via the "
                                                                         "'Source & Epoch' panel to "
                                                                         "compare them side-by-side.")])),
        # Container that the compute callback fills with either a simple panel
        # (when no target sources are specified) or a tabbed view (one tab per source).
        html.Div(id='outputs-container', className='m-0 p-0')
    ])
