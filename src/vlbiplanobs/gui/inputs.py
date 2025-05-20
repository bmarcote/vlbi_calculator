import functools
from typing import Optional, Sequence
from datetime import datetime as dt
from dash import html, dcc
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from vlbiplanobs import freqsetups as fs
from vlbiplanobs import stations
from vlbiplanobs import observation


def modal_welcome() -> html.Div:
    """Welcoming window announcing the new version
    """
    return html.Div([
        dcc.Store(id='welcome-modal-shown', storage_type='local'),
        dbc.Modal([
            dbc.ModalBody([
                    html.H4("Welcome to the new EVN Observation Planner (PlanObs)!"),
                    html.H6("The version 3.0 (beta) is here"),
                    html.P("This version contains a major upgrade of PlanObs and a big re-design."
                           "The main new features are:"),
                    html.Ul([
                        html.Li("Real antenna mount limits are now taken into account to know if a "
                                "source can be observed (e.g. hour angle limits)."),
                        html.Li("The sensitivity calculations take into account the limited bandwidth "
                                "that some antennas may have (e.g. eMERLIN stations observing within an "
                                "EVN observation)."),
                        html.Li("Warnings if the Sun is too close to your source during the observation."),
                        html.Li("Fully-featured command-line program (CLI) that allows you to plan "
                                "observations locally, and with multiple sources or from catalogs ("
                                "this will be added to the online PlanObs in a later release)."),
                        html.Li("The speed of PlanObs has increased significantly!"),
                        html.Li("And yet, several new features are yet to come in the next releases.")]),
                    html.H3(html.Em("Enjoy!"))], style={'background-color': 'rgb(212, 212, 216)'}),
            dbc.ModalFooter(dbc.Button("Close", id="close-modal", className="ml-auto btn btn-secondary"))],
            id="welcome-modal", is_open=False)])


def modal_general_info() -> html.Div:
    """Returns the modal window that shows all relevant information concerning PlanObs.
    """
    return html.Div(className='sidebar',
                    children=[
                        html.H5("EVN Observation Planner"),
                        html.P(["The ", html.B('EVN Observation Planner'), " (", html.Em('PlanObs'),
                                ") is a tool that allows users to plan very long baseline "
                                "interferometry (VLBI) observations and verify the "
                                "expectations of such observations."]),
                        html.P(["PlanObs is mainly developed and maintained by ",
                                html.A("Benito Marcote", href="https://bmarcote.github.io/",
                                       target="_blank"),
                                " at the ", html.A("Joint Institute for VLBI ERIC (JIVE)",
                                                   href="https://www.jive.eu", target="_blank"),
                                " to primarily support observations with the ",
                                html.A("European VLBI Network (EVN)", href="https://www.evlbi.org",
                                       target="_blank"),
                                ". PlanObs's code is publicly available at ",
                                html.A("GitHub", href="https://github.com/bmarcote/vlbi_calculator",
                                       target="_blank"),
                                " under the GNU GPLv3.0+ license."
                                ]),

                        html.Hr(className='horizontal mb-1 d-xl-block d-none dark'),
                        html.H6("Command-Line Interface"),
                        html.P([html.B("PlanObs"), " can be used locally in your computer as a "
                                "command-line program with different capabilities and options that "
                                "extend the functionality of the web interface. "]),
                        html.P(["You can install PlanObs directly from ", html.B("pip"), " as ", html.Br(),
                                html.Span("> python3 -m pip install vlbiplanobs",
                                          className='text-center',
                                          style={'color': "#004990"})]),
                        html.Div(className='text-center col-12', children=[
                            dbc.Button("Read the Documentation",
                                       href="https://github.com/bmarcote/vlbi_calculator/blob/master"
                                            "/README.md", target="_blank", external_link=True,
                                       className='btn btn-secondary', outline=True)]),

                        html.Hr(className='horizontal mb-1 d-xl-block d-none dark'),
                        html.H6("Issues or Feature Requests"),
                        html.P(["PlanObs is currently undergoing a major update. It is thus not "
                                "unlikely that you may encounter some small bugs or missing features. "
                                "In such cases, or either if you just need further information on the "
                                "current capabilities, or you would like to see new features, please "
                                "do not heasitate to contact us! "]),
                        html.Div(className='text-center',
                                 children=dbc.Button(
                                    "Open an Issue",
                                    href="https://github.com/bmarcote/vlbi_calculator/issues",
                                    target="_blank", external_link=True,
                                    className='btn btn-secondary', outline=True))
                    ])


def card(children: Optional[list | html.Div] = None, className: Optional[str] = None,
         classNameDiv: Optional[str] = None, style: dict = {}) -> html.Div:
    return html.Div(className='m-0 p-2' if classNameDiv is None else classNameDiv,
                    children=html.Div(className='card m-0 p-0' if className is None else 'card ' + className,
                                      style={'margin': '10px'} if not style else style,
                                      children=html.Div(className='card-body', children=children)))


def antenna_card_hover(app, target, ant: stations.Station, show_wavelengths: bool = False) -> html.Div:
    """Creates a Card showing all the relevant information for a given antenna, which
    will be displayed as a pop up window when hovering the name of a given antenna.
    """
    return html.Div(dmc.HoverCard(shadow="md", withArrow=True, width=300,
                                  openDelay=300, position='right', classNames='col-12 p-0 m-0',
                    children=[dmc.HoverCardTarget(target),
                              dmc.HoverCardDropdown([antenna_card(app, ant, show_wavelengths)],
                                                    className='col-12 m-0 p-0')]),
                    style={'display': 'inline-flex', 'flex-wrap': 'wrap',
                           'grid-template-columns': 'repeat(auto-fit, minmax(10rem, 1fr))'})


def parse_str_list(a_list: Sequence[str]) -> str:
    """Returns a string with the list elements as a comma-separated string list.
    For the last element, it adds an 'and' instead of the comma.
    """
    return ', '.join(a_list)[::-1].replace(',', ' dna ,' if len(a_list) > 3 else 'dna ', 1)[::-1]


def print_table_bands_sefds(ant: stations.Station, show_wavelengths: bool = False) -> html.Div:
    bands_str = [fs.bands[band].split('or')[0].replace('cm', ' cm').strip() if show_wavelengths else
                 fs.bands[band].split('or')[1].replace('GHz', ' GHz') for band in ant.bands]

    return html.Div([html.Div(className='col-6 px-0 m-0',
                              children=[
                                html.Label(children=b_str,
                                           className='text-xs text-primary'),
                                html.Label(f"({ant.sefds[band].value:g} "
                                           f"{ant.sefds[band].unit.to_string('unicode')})",
                                           className='text-xs text-secondary')]
                              ) for b_str, band in zip(bands_str, ant.bands)], className='row')


def antenna_card(app, ant: stations.Station, show_wavelengths: bool = True) -> html.Div:
    return html.Div(dmc.Card(children=[html.Div(
                        dmc.CardSection(
                            dmc.Image(src=app.get_asset_url(
                                f"ant-{ant.name.replace(' ', '_').lower()}.jpg"), h='300px', alt=ant.name),
                        ), className='col-12 p-0 m-0'),
                        html.Div([
                            dmc.Group([
                                dmc.Text(f"{ant.name} ({ant.codename})", className='text-bolder'),
                                dmc.Badge(ant.diameter, color='#9DB7C4', style={'text-transform': 'none'}),
                            ], justify='space-between', mt='md', pt='1rem', mb='0'),
                            dmc.Text(ant.fullname, mb='0', c='#004990') if ant.fullname != ant.name
                            else None,
                            dmc.Text(ant.country, c='dimmed', mt='0', mb='1rem'),
                            dmc.Text(f"Default antenna in {parse_str_list(ant.networks)}.", mb='1rem',
                                     size='sm') if len(ant.networks) > 0 else None,
                            dmc.Text("No longer operational.", mb='1rem', c='#a01d26', size='sm')
                            if ant.decommissioned else None,
                            dmc.Text("Can observe at the following bands (System Equivalent Flux "
                                     "Density, SEFD, values in brackets):", size='sm',
                                     id=f"ant-{ant.codename}-band-spec"),
                            html.Div(className='container col-12 row mx-0 px-0',
                                     id=f"badge-band-ant-{ant.codename}",
                                     children=print_table_bands_sefds(ant, show_wavelengths))
                        ], className='col-12 px-2 pb-2 m-0')
                    ], withBorder=False), className='col-12 p-0 m-0', style={'width': '300px'})


def compute_button() -> html.Div:
    """Returns the button to compute the observation
    """
    return html.Div(
                html.Div([
                    html.Button('CALCULATE',
                                id='compute-observation',
                                className='btn btn-evn text-bolder btn-lg mx-auto w-75 m-4 p-2',
                                style={'position': 'sticky', 'top': '20px'}),
                    dbc.Spinner(id='loading', color='#a01d26',
                                children=html.Div(id='loading-div'))
                ], className='d-flex align-items-center justify-content-center', style={'gap': '5px'}),
        className='col-6', style={'position': 'relative'})


@functools.cache
def pick_band_labels(show_wavelengths: bool) -> list:
    if show_wavelengths:
        return [{'label': 'λ\n(cm)', 'style': {'font-weight': 'bold', 'color': '#004990',
                                               'white-space': 'pre'}}] + \
               [b.split('or')[0].replace('cm', '').strip() for b in fs.bands.values()]

    return [{'label': 'ν\n(GHz)', 'style': {'font-weight': 'bold', 'color': '#004990',
                                            'white-space': 'pre'}}] + \
           [b.split('or')[1].replace('GHz', '').strip() for b in fs.bands.values()]


@functools.cache
def band_from_index(ind: int) -> str:
    """Returns the band selected by the slider, given the index of the slider.
    If no band is selected, returns None.
    """
    if ind == 0:
        raise ValueError("No band selected")

    return list(fs.bands.keys())[ind - 1]


def pick_band(bands: dict[str, str]) -> html.Div:
    """Returns the card allowing the user to pick up the band.
    """
    labels = pick_band_labels(False)
    return html.Div(className='row col-12',
                    children=[html.Div(className='col-12', children=[
                       html.Div(className='row align-items-bottom', children=[
                                html.Div(className='col-8', children=[
                                  html.H4("Observing band", className='text-dark font-weight-bold mb-0'),
                                ]),
                                html.Div(className='col-4 text-right', children=[
                                   dbc.Switch(label='Show wavelengths', value=False, id='switch-band-label',
                                              persistence=True),
                                ])]),
                       html.Br(),
                       html.Div(dcc.Slider(min=0, max=len(labels), step=1, value=0,
                                           marks={i: label for i, label in enumerate(labels)},
                                included=False, id='band-slider', persistence=True))])])


def network_band_labels(network: str, show_wavelengths: bool = False) -> str:
    if show_wavelengths:
        label_bands = [b.replace('cm', '').strip() for b in observation._NETWORKS[network].observing_bands]
    else:
        label_bands = [fs.bands[b].split('or')[1].replace('GHz', '').strip()
                       for b in observation._NETWORKS[network].observing_bands]

    if len(label_bands) == 0:
        return 'N/A'
    elif len(label_bands) == 1:
        return label_bands[0]
    else:
        return ', '.join(label_bands[:-1]) + (', and ' if len(label_bands) > 3 else ' and ') + \
            label_bands[-1] + (' cm.' if show_wavelengths else ' GHz.')


def network_entry(app, network: str) -> html.Div:
    """Returns the DIV that shows a given VLBI network to be selected
    """

    has_full_name = stations.Stations.get_network_full_name(network) != network
    return html.Div(dbc.Card([
        dbc.CardImg(src=app.get_asset_url(f"network-{network.lower()}.png"),
                    top=False, className='full-background',
                    style={'opacity': 1.0, 'object-fit': 'cover', 'height': '10rem'}),
        dbc.CardImgOverlay(dbc.CardBody([
            html.Div(className='d-flex align-items-center', children=[
                html.H4(network, className='text-bold text-white mb-0 mt-0 flex-grow-1'),
                dbc.Switch(label='', value=False, id=f"network-{network}",
                           className='form-check ml-auto', persistence=True)]),
            html.H6(stations.Stations.get_network_full_name(network),
                    className='mt-0 mb-0 text-white',
                    style={'white-space': 'nowrap', 'text-overflow': 'ellipsis',
                           'overflow': 'hidden'}) if has_full_name else html.Br(),
            html.Label(f"Observes at {network_band_labels(network)}",
                       className='card-text text-white',
                       id=f"network-{network}-label-band",
                       style={'line-height': '1.2'})],
                       className='card-body text-start p-0 pt-0 w-100'))],
                              className='m-2 card-background',
                              style={'min-width': '11rem', 'height': '8rem',
                                     'overflow': 'hidden', 'opacity': 1.0},
                              id=f"network-{network}-card"), className='col-4')


def networks(app) -> html.Div:
    return html.Div([html.H4("Select default VLBI Network(s)",
                             className='text-dark font-weight-bold mb-0 pl-2 ml-4'),
                     html.Div([network_entry(app, network_name) for network_name in observation._NETWORKS],
                              className='d-flex flex-wrap')])


def antenna_list(app, show_wavelengths: bool = False) -> html.Div:
    """Returns the DIV that shows all the antennas that can be selected, and allows searching through them
    """
    return html.Div(dbc.Accordion(start_collapsed=True, className='accordion-contents-open-ids', flush=True,
                                  persistence=True, id='accordion-state',
                                  children=dbc.AccordionItem(
                                      title=html.H4(['Manual Selection of Antennas   ',
                                                     html.I(className='fa fa-solid fa-angle-down',
                                                            id='accordion-ant')],
                                                    className='text-dark font-weight-bold mb-0 '
                                                    'accordion-header'),
                                      children=[dmc.Group(
                                          dmc.ChipGroup(value=[], id='switches-antennas', persistence=True,
                                                        multiple=True, deselectable=True, children=[
                                              antenna_card_hover(app,
                                                                 dmc.Chip(s.name, value=s.codename,
                                                                          id=f"chip-{s.codename}",
                                                                          color='#004990',
                                                                          persistence=True,
                                                                          styles={'display': 'grid',
                                                                                  'grid-template-columns':
                                                                                  'repeat(auto-fit, '
                                                                                  'minmax(10rem, 1fr))'}),
                                                                 s, show_wavelengths)
                                              for s in observation._STATIONS]),
                                          className='container mb-2 flex',
                                          style={'display': 'inline-flex', 'gap': '5px',
                                                 'flex-wrap': 'wrap'})])))


def source_selection() -> html.Div:
    return html.Div([html.H4("Source To Observe", className='text-dark font-weight-bold mb-1'),
                     dbc.Switch(label='Specify a name or coordinates', value=False,
                                id='switch-specify-source', persistence=True),
                     html.Div(id='source-selection-div', hidden=True,
                              children=html.Div(className='row', children=[
                                dcc.Input(id='source-input', value=None, type='text',
                                          className='form-control', placeholder="hh:mm:ss dd:mm:ss",
                                          persistence=True, debounce=True),
                                html.Small(id='error_source', className='form-text text-muted')])),
                     html.Br(),
                     html.Div(className='row form-group', children=[
                         html.Label('Percentage of observing time spent on target', htmlFor='onsourcetime'),
                         dcc.Slider(id='onsourcetime', min=20, max=100, step=10, value=70,
                                    marks={i: str(i) for i in range(20, 101, 10)},
                                    persistence=True)])])


def epoch_selection() -> html.Div:
    return html.Div([html.H4("Define Observing Epoch", className='text-dark font-weight-bold mb-1'),
                     dbc.Switch(label='Specify epoch', value=False,
                                id='switch-specify-epoch', persistence=True),
                     html.Div(id='epoch-selection-div', hidden=True, children=[
                        html.Div(className='row', children=[
                            html.Label('Start of observation (UTC)', htmlFor='starttime'),
                            html.Div(className='row mx-0', children=[
                                html.Div(className='col-12 px-0 mx-0', children=[
                                    dmc.DateTimePicker(id='starttime', className='form-picker',
                                                       value=None, minDate=dt(1900, 1, 1),
                                                       maxDate=dt(2100, 1, 1),
                                                       valueFormat='DD-MM-YYYY HH:mm',
                                                       placeholder='Start time', withSeconds=False,
                                                       persistence=True, clearable=True,
                                                       style={'width': '100%'})])]),
                            html.Div(className='row', children=[
                                html.Small(id='error_starttime', style={'color': 'red'},
                                           className='form-text text-muted')])])]),
                    html.Div(className='row form-group', children=[
                        html.Label('Duration of the observation (in hours)', htmlFor='duration'),
                        dcc.Input(id='duration', value=None, type='number', className='form-control',
                                  placeholder="Duration of the observation (in hours)",
                                  persistence=True, debounce=True, inputMode='numeric', min=0.5, max=50),
                        html.Small(id='error_duration', className='form-text text-muted')])])


def correlations() -> html.Div:
    """Creates the inputs to specify if the user wants continuum correlation and/or spectral line.
    """
    return html.Div([html.H4("Observation Setup", className='text-dark font-weight-bold mb-1'),
                     html.Div(className='col-12', children=[
                        html.Div(className='row d-flex align-items-bottom', children=[
                            html.Div(className='col-6', children=[
                                dbc.Switch(label='Real-time (e-EVN observation)', value=False,
                                           id='switch-specify-e-evn', persistence=True),
                                dbc.Switch(label='Optimize for spectral line', value=False,
                                           id='switch-specify-continuum', persistence=True),
                                dcc.Dropdown(id='datarate', placeholder="Maximum observing data rate...",
                                             options=tuple({'label': drl, 'value': str(dr)}
                                                           for dr, drl in fs.data_rates.items()),  # type: ignore
                                             value=2048, persistence=True, clearable=False),
                                html.Label(id='bandwidth-label', style={'color': '#999999'}, children='',
                                           htmlFor="datarate")]),
                            html.Div(className='col-6', children=[
                                html.Div(className='form-group', children=[
                                    dcc.Dropdown(id='subbands', placeholder="Select no. subbands...",
                                                 options=tuple({'label': sbl, 'value': sb}  # type: ignore
                                                               for sb, sbl in fs.subbands.items())[::-1],
                                                 value=8, persistence=True, clearable=False),
                                    dcc.Dropdown(id='channels',
                                                 placeholder="Select no. channels per subband...",
                                                 options=tuple({'label': chl, 'value': ch}  # type: ignore
                                                               for ch, chl in fs.channels.items())[::-1],
                                                 value=64, persistence=True, clearable=False),
                                    dcc.Dropdown(id='pols', placeholder="Polarizations...",
                                                 options=tuple({'label': pl, 'value': p}  # type: ignore
                                                               for p, pl in fs.polarizations.items()),
                                                 value=4, persistence=True, clearable=False),
                                    dcc.Dropdown(id='inttime', placeholder="Integration time...",
                                                 options=tuple({'label': itl, 'value': it}  # type: ignore
                                                               for it, itl in fs.inttimes.items()),
                                                 value=2, persistence=True, clearable=False)])])])])])
