from dash import html, dcc
import dash_bootstrap_components as dbc
from vlbiplanobs import freqsetups as fs
from vlbiplanobs.gui import inputs, outputs


"""Module that assembles the top-level page layout blocks (banner, input column, compute
buttons, output column) from the individual Dash components defined in `gui.inputs` and
`gui.outputs`. `gui.main` wires these blocks together into the final `app.layout`.
"""


def top_banner(app) -> html.Div:
    """Builds the top banner card with the EVN and JIVE logos and the page title.

    Parameters
    ----------
    app : dash.Dash
        The Dash app instance, used to resolve asset URLs (`app.get_asset_url`) for the logos.

    Returns
    -------
    html.Div
        The banner as a `html.Div` (wrapped in an `inputs.card`).
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
    """Builds the left-hand column of input cards: band, duration, source/epoch, networks,
    antenna list, and correlation setup.

    Parameters
    ----------
    app : dash.Dash
        The Dash app instance, forwarded to the `inputs.networks` and `inputs.antenna_list`
        builders that need it to resolve asset URLs.

    Returns
    -------
    html.Div
        The input column as a `html.Div` containing one `inputs.card` per input section.
    """
    return html.Div(children=[
            inputs.card(inputs.pick_band(fs.bands)),
            inputs.card(inputs.duration()),
            inputs.card(inputs.source_and_epoch_selection()),
            inputs.card(inputs.networks(app)),
            inputs.card(inputs.antenna_list(app)),
            inputs.card(inputs.correlations())])


def compute_buttons(app) -> html.Div:
    """Builds the compute/download button row for the non-real-time (button-triggered) mode.

    Parameters
    ----------
    app : dash.Dash
        The Dash app instance. Currently unused by this function, kept for signature
        consistency with the other layout builders.

    Returns
    -------
    html.Div
        The button row as a `html.Div` containing the compute button, the download button, and
        the `dcc.Download` component.
    """
    return html.Div([html.Div(className='m-0 p-0', children=[
        html.Div(className='row d-flex m-0 p-0', children=[
            inputs.compute_button(),
            outputs.download_button_div(),
            dcc.Download(id="download-data")])]),
            # html.Div(className='m-0 p-0', children=[
            #     dcc.Location(id='url', refresh=False),
            #     inputs.export_button(),
            #     Alert("Nothing to see here", id='export-alert',
            #           is_open=False, color='success',
            #           duration=5000, dismissable=True)
            # ])
    ])


def compute_buttons_realtime(app) -> html.Div:
    """Layout for real-time mode: just the loading spinner.

    The PDF download button is rendered per-target inside each output tab (and inside
    the simple panel when no target source is specified), so a global button is no
    longer needed here.

    Parameters
    ----------
    app : dash.Dash
        The Dash app instance. Currently unused by this function, kept for signature
        consistency with the other layout builders.

    Returns
    -------
    html.Div
        The spinner row as a `html.Div`.
    """
    return html.Div(className='m-0 p-0', children=[
        html.Div(className='row d-flex m-0 p-0 justify-content-center', children=[
            html.Div(className='col-12 text-center', children=[
                dbc.Spinner(id='loading', color='#a01d26',
                            children=html.Div(id='loading-div'))])])])


def outputs_column(app) -> html.Div:
    """Builds the right-hand output column: the shared top-level user message plus the empty
    `outputs-container` that the compute callback fills in with the results.

    Parameters
    ----------
    app : dash.Dash
        The Dash app instance. Currently unused by this function, kept for signature
        consistency with the other layout builders.

    Returns
    -------
    html.Div
        The output column as a `html.Div` with the `user-message` and `outputs-container` children.
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
