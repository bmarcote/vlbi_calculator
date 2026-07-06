import functools
import importlib
from typing import Optional, Sequence
from datetime import datetime as dt
from dash import html, dcc
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from vlbiplanobs import freqsetups as fs
from vlbiplanobs import stations
from vlbiplanobs import observation

"""Module that builds the Dash layout components (inputs) shown in the PlanObs GUI:
modals, antenna/network selection cards, band/duration/correlation pickers, and the
target-source modal. These are pure layout-construction functions (no callbacks);
the interactive behavior is wired up in `vlbiplanobs.gui.callbacks`.
"""


def modal_welcome() -> html.Div:
    """Welcoming window announcing the new version.
    """
    return html.Div([
        dcc.Store(id='welcome-modal-shown', storage_type='local'),
        dbc.Modal([
            dbc.ModalBody([
                    html.H4("Welcome to the new EVN Observation Planner (PlanObs) version!"),
                    html.H6(f"The version {importlib.metadata.version('vlbiplanobs')} is here"),
                    html.P("We have been fixing bugs that came with the new major version "
                           "and you may have noticed in the last month."),
                    html.Ul([
                        html.Li("The PDF summary should now give you the correct information."),
                        html.Li("Small bug fixed in the elevation plots."),
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

    Returns
    -------
    html.Div
        The modal window content.
    """
    return html.Div(className='sidebar',
                    children=[
                        html.H5("EVN Observation Planner", className='mb-0 pb-0'),
                        html.H6(f"(version {importlib.metadata.version('vlbiplanobs')})", className='mt-0 pt-0'),
                        # Theme toggle button
                        html.Button(
                            id='theme-toggle-btn',
                            className='theme-toggle-btn',
                            children=[
                                html.I(className='fas fa-sun theme-icon-light'),
                                html.I(className='fas fa-moon theme-icon-dark'),
                                html.Span('Toggle Dark Mode', id='theme-toggle-text')
                            ]
                        ),
                        html.Hr(className='horizontal mb-1 d-xl-block d-none dark'),
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
         classNameDiv: Optional[str] = None, style: Optional[dict] = None) -> html.Div:
    """Returns a generic Bootstrap-style card wrapper (outer DIV, card DIV, card-body DIV)
    around the given children.

    Parameters
    ----------
    children : list or html.Div, optional
        Content to place inside the card body. None by default.
    className : str, optional
        Extra CSS class(es) appended to 'card' for the card DIV. None by default.
    classNameDiv : str, optional
        CSS class for the outer wrapper DIV. Defaults to 'm-0 p-2' if not given.
    style : dict, optional
        Inline style for the card DIV. Defaults to {'margin': '10px'} if not given.

    Returns
    -------
    html.Div
        The card wrapper DIV.
    """
    return html.Div(className='m-0 p-2' if classNameDiv is None else classNameDiv,
                    children=html.Div(className='card m-0 p-0' if className is None else 'card ' + className,
                                      style={'margin': '10px'} if style is None else style,
                                      children=html.Div(className='card-body', children=children)))


def antenna_card_hover(app, target, ant: stations.Station, show_wavelengths: bool = False) -> html.Div:
    """Creates a Card showing all the relevant information for a given antenna, which
    will be displayed as a pop up window when hovering the name of a given antenna.

    Parameters
    ----------
    app : Dash app instance
        Used to resolve the antenna image asset URL.
    target : Dash component
        The component that triggers the hover card (e.g. an antenna chip).
    ant : stations.Station
        The antenna to describe in the card.
    show_wavelengths : bool, optional
        If True, shows the observing bands in wavelength (cm); otherwise in frequency (GHz).
        False by default.

    Returns
    -------
    html.Div
        The hover-card wrapper DIV.
    """
    return html.Div(dmc.HoverCard(shadow="lg", radius="lg",
                                  openDelay=700, position='right',
                    children=[dmc.HoverCardTarget(target),
                              dmc.HoverCardDropdown(antenna_card(app, ant, show_wavelengths),
                                                    className='m-0 p-0')]),
                    style={'display': 'inline-flex', 'flex-wrap': 'wrap',
                           'grid-template-columns': 'repeat(auto-fit, minmax(10rem, 1fr))'})


def parse_str_list(a_list: Sequence[str]) -> str:
    """Returns a string with the list elements as a comma-separated string list.
    For the last element, it adds an 'and' instead of the comma.

    Parameters
    ----------
    a_list : Sequence[str]
        The list of strings to join.

    Returns
    -------
    str
        A single string with all elements joined, e.g. "A, B, and C".
    """
    return ', '.join(a_list)[::-1].replace(',', ' dna ,' if len(a_list) > 3 else 'dna ', 1)[::-1]


def print_table_bands_sefds(ant: stations.Station, show_wavelengths: bool = False) -> html.Div:
    """Returns a two-column DIV listing, for each band the antenna can observe, the band label
    (wavelength or frequency) alongside its SEFD value.

    Parameters
    ----------
    ant : stations.Station
        The antenna whose observable bands and SEFDs are listed.
    show_wavelengths : bool, optional
        If True, labels bands in wavelength (cm); otherwise in frequency (GHz). False by default.

    Returns
    -------
    html.Div
        The two-column DIV with band labels and SEFD values.
    """
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
    """Returns the Card content describing an antenna: picture, name, diameter, country,
    default network membership, decommissioned status, and observable bands with SEFDs.

    Parameters
    ----------
    app : Dash app instance
        Used to resolve the antenna image asset URL.
    ant : stations.Station
        The antenna to describe.
    show_wavelengths : bool, optional
        If True, shows the observing bands in wavelength (cm); otherwise in frequency (GHz).
        True by default.

    Returns
    -------
    html.Div
        The Card content DIV.
    """
    return html.Div(dmc.Card(children=[html.Div(
                        dmc.CardSection(
                            dmc.Image(src=app.get_asset_url(
                                f"ant-{ant.name.replace(' ', '_').lower()}.jpg"), h='200px', alt=ant.name),
                        ), className='col-12 p-0 m-0'),
                        html.Div([
                            dmc.Group([
                                dmc.Text(f"{ant.name} ({ant.codename})", className='text-bolder', size='sm'),
                                dmc.Badge(ant.diameter, color='#9DB7C4', style={'text-transform': 'none'}),
                            ], justify='space-between', mt='md', pt='0.5rem', mb='0'),
                            dmc.Text(ant.fullname, mb='0', c='#004990', size='sm') if ant.fullname != ant.name
                            else None,
                            dmc.Text(ant.country, c='dimmed', mt='0', mb='1rem'),
                            dmc.Text(f"Default antenna in {parse_str_list(ant.networks)}.", mb='1rem',
                                     size='sm') if len(ant.networks) > 0 else None,
                            dmc.Text("No longer operational.", mb='1rem', c='#a01d26', size='sm')
                            if ant.decommissioned else None,
                            dmc.Text("Can observe at the following bands (System Equivalent Flux "
                                     "Density, SEFD, values in brackets):", size='sm'),
                            html.Div(className='container col-12 row mx-0 px-0 text-xs',
                                     id={'type': 'badge-band-ant', 'index': ant.codename},
                                     children=print_table_bands_sefds(ant, show_wavelengths))
                        ], className='col-12 px-1 pb-2 m-0')
                    ], withBorder=False), className='col-12 p-0 m-0 text-sm', style={'width': '200px'})


def compute_button() -> html.Div:
    """Returns the button to compute the observation.
    """
    return html.Div(
                html.Div([
                    html.Button('CALCULATE',
                                id='compute-observation',
                                className='btn btn-evn text-bolder btn-lg mx-auto w-75 m-4 p-2',
                                style={'position': 'sticky', 'top': '20px', 'z-index': '1000'}),
                    dbc.Spinner(id='loading', color='#a01d26',
                                children=html.Div(id='loading-div'))
                ], className='d-flex align-items-center justify-content-center', style={'gap': '5px'}),
        className='col-6', style={'position': 'relative'})

def export_button() -> html.Div:
    """Returns the button to export the observation to Polaris.
    """
    return html.Div(
                html.Div([
                    html.Button('Export to Polaris',
                                id='export-state-of-the-system',
                                className='btn btn-evn text-bolder btn-lg mx-auto w-75 m-4 p-2',
                                style={'position': 'sticky', 'top': '20px', 'z-index': '1000'}),
                ], className='d-flex align-items-center justify-content-center', style={'gap': '5px'}),
        className='col-6', style={'position': 'relative'})


@functools.cache
def pick_band_labels(show_wavelengths: bool) -> list:
    """Returns the list of slider marks/labels for the band-selection slider, prefixed with
    a bold axis-name label ('λ (cm)' or 'ν (GHz)').

    Cached since the result only depends on `show_wavelengths`.

    Parameters
    ----------
    show_wavelengths : bool
        If True, labels the bands by wavelength (cm); otherwise by frequency (GHz).

    Returns
    -------
    list
        The axis label dict followed by one label string per band, in `freqsetups.bands` order.
    """
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

    Parameters
    ----------
    ind : int
        The index of the slider.

    Returns
    -------
    str
        The band key selected by the slider.

    Raises
    ------
    ValueError
        If no band is selected (index is 0).
    """
    if ind == 0:
        raise ValueError("No band selected")

    return list(fs.bands.keys())[ind - 1]


def pick_band(bands: dict[str, str]) -> html.Div:
    """Returns the card allowing the user to pick up the band.

    Parameters
    ----------
    bands : dict[str, str]
        Currently unused by the function body (labels are always built from
        `freqsetups.bands` via `pick_band_labels`); kept for API compatibility.

    Returns
    -------
    html.Div
        The band-picker card DIV.
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
                       html.Div(dcc.Slider(min=0, max=len(labels)-1, step=1, value=0,
                                           marks={i: label for i, label in enumerate(labels)},
                                included=False, id='band-slider', persistence=True,
                                allow_direct_input=False))])])


def network_band_labels(network: str, show_wavelengths: bool = False) -> str:
    """Returns a human-readable, comma-separated string of the bands a network can observe.

    Parameters
    ----------
    network : str
        Key of the network in `observation._NETWORKS`.
    show_wavelengths : bool, optional
        If True, lists bands in wavelength (cm); otherwise in frequency (GHz). False by default.

    Returns
    -------
    str
        'N/A' if the network has no observing bands, the single band label if there is
        only one, or a comma-separated list ending in 'and <last band> cm.'/'GHz.' otherwise.
    """
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
    """Returns the DIV that shows a given VLBI network to be selected.

    Parameters
    ----------
    app : Dash app instance
        Used to resolve the network background image asset URL.
    network : str
        Key of the network in `observation._NETWORKS`.

    Returns
    -------
    html.Div
        The network card DIV.
    """

    # has_full_name = stations.Stations.get_network_full_name(network) != network
    return html.Div(dbc.Card([
        dbc.CardImg(src=app.get_asset_url(f"network-{network.lower()}.png"),
                    top=False, className='full-background',
                    style={'opacity': 0.5, 'object-fit': 'cover', 'height': '7rem'}),
        dbc.CardImgOverlay(dbc.CardBody([
            html.Div(className='d-flex align-items-center', children=[
                html.H4(network, className='text-bold text-white px-1 mb-0 mt-0 flex-grow-1',
                        title=stations.Stations.get_network_full_name(network)),
                dbc.Switch(label='', value=False, id={'type': 'network-switch', 'index': network},
                           className='form-check ml-auto', persistence=True)]),
            # html.Br(),
            # html.H6(stations.Stations.get_network_full_name(network),
            #         className='mt-0 mb-0 text-white',
            #         style={'white-space': 'nowrap', 'text-overflow': 'ellipsis',
            #                'overflow': 'hidden'}) if has_full_name else html.Br(),
            html.Label(network_band_labels(network),
                       className='card-text text-white mt-auto',
                       id={'type': 'network-label-band', 'index': network},
                       style={'line-height': '1.2'})],
                       className='card-body text-start p-0 pt-0 w-100 h-100 d-flex flex-column'), className='p-2')],
                              className='m-2 card card-background',
                              style={'min-width': '11rem', 'height': '7rem',
                                     'overflow': 'hidden', 'opacity': 1.0},
                              id={'type': 'network-card', 'index': network}), className='col-4')


def networks(app) -> html.Div:
    """Returns the DIV listing all default VLBI networks, each rendered as a selectable card
    (via `network_entry`).

    Parameters
    ----------
    app : Dash app instance
        Passed through to `network_entry` to resolve image asset URLs.

    Returns
    -------
    html.Div
        The DIV listing all default VLBI network cards.
    """
    return html.Div([html.H4("Select default VLBI Network(s)",
                             className='text-dark font-weight-bold mb-0 pl-2 ml-4'),
                     html.Div([network_entry(app, network_name) for network_name in observation._NETWORKS],
                              className='d-flex flex-wrap')])


def station_groups() -> dict[str, list[stations.Station]]:
    """Returns an ordered dict mapping group_name -> [Station, ...] for all grouped stations.
    Stations without a group are not included. Order follows the catalog order.

    Returns
    -------
    dict[str, list[stations.Station]]
        Mapping of group name to the list of stations in that group.
    """
    groups: dict[str, list[stations.Station]] = {}
    for s in observation._STATIONS:
        if s.group is not None:
            groups.setdefault(s.group, []).append(s)
    return groups


def _grouped_chip_component(app, group_name: str, group_stations: list[stations.Station],
                             show_wavelengths: bool = False) -> html.Div:
    """Renders a single grouped-antenna chip: a toggle button (on/off) plus a small dropdown
    arrow that lets the user pick which configuration is active. Only one config per group can
    be in the selected set at a time. Each dropdown item shows a HoverCard tooltip with the
    antenna info card.

    Parameters
    ----------
    app : Dash app instance
        Used to resolve antenna image asset URLs.
    group_name : str
        Group identifier (e.g. 'vla').
    group_stations : list[stations.Station]
        Stations in this group, in catalog order.
    show_wavelengths : bool, optional
        Passed through to `antenna_card`. False by default.

    Returns
    -------
    html.Div
        The grouped-antenna chip wrapper DIV.
    """
    default_codename = group_stations[0].codename

    stores = [
        dcc.Store(id={'type': 'group-active-codename', 'index': group_name}, data=default_codename,
                  storage_type='session'),
        dcc.Store(id={'type': 'group-is-selected', 'index': group_name}, data=False,
                  storage_type='session'),
    ]

    menu_items = []
    for i, s in enumerate(group_stations):
        item_content = dmc.MenuItem(
            s.name,
            id={'type': 'group-menu-item', 'index': f"{group_name}__{s.codename}"},
            style={'font-size': '0.85rem'},
            className='group-menu-item group-menu-item-active' if i == 0 else 'group-menu-item',
        )
        menu_items.append(
            dmc.HoverCard(shadow="lg", radius="lg", openDelay=700, position='right',
                          children=[dmc.HoverCardTarget(item_content),
                                    dmc.HoverCardDropdown(antenna_card(app, s, show_wavelengths),
                                                          className='m-0 p-0')])
        )

    dropdown_trigger = dmc.Menu(
        id={'type': 'group-menu', 'index': group_name},
        position='bottom-start',
        children=[
            dmc.MenuTarget(
                html.Button(
                    '▼',
                    id={'type': 'group-dropdown-btn', 'index': group_name},
                    className='btn-group-config-arrow',
                    title='Switch configuration',
                )
            ),
            dmc.MenuDropdown(menu_items),
        ]
    )

    # The toggle button (on/off) shows the active configuration name (e.g. 'VLA 1').
    # Clicking it toggles selection; the label is kept in sync by a callback.
    toggle_btn = html.Button(
        id={'type': 'group-toggle-btn', 'index': group_name},
        children=group_stations[0].name,
        className='btn-group-chip-toggle btn-group-chip-off',
        title=f"Toggle {group_name.upper()} antenna",
    )

    wrapper = html.Div(
        stores + [
            html.Div(
                [toggle_btn, dropdown_trigger],
                className='group-chip-inner',
            )
        ],
        id={'type': 'group-chip-wrapper', 'index': group_name},
        style={'display': 'inline-flex', 'align-items': 'center', 'align-self': 'center'}
    )
    return wrapper


def antenna_list(app, show_wavelengths: bool = False) -> html.Div:
    """Returns the DIV that shows all the antennas that can be selected, and allows searching
    through them.

    Ungrouped antennas are rendered as individual chips inside the `switches-antennas`
    `dmc.ChipGroup`; antennas that belong to a group (see `station_groups`) are instead
    rendered once per group as a combined toggle+dropdown chip via `_grouped_chip_component`.

    Parameters
    ----------
    app : Dash app instance
        Passed through to build antenna hover cards.
    show_wavelengths : bool, optional
        If True, shows the observing bands in wavelength (cm); otherwise in frequency (GHz).
        False by default.

    Returns
    -------
    html.Div
        The DIV listing all selectable antenna chips.
    """
    groups = station_groups()
    grouped_codenames = {s.codename for slist in groups.values() for s in slist}

    grouped_components = [
        _grouped_chip_component(app, gname, gstations, show_wavelengths)
        for gname, gstations in groups.items()
    ]

    return html.Div([html.H4("Manual Selection of Antennas   ",
                             className='text-dark font-weight-bold mb-2 pl-2 ml-4'),
                     dmc.Group(
                         [
                             dmc.ChipGroup(value=[], id='switches-antennas', persistence=True,
                                           multiple=True, deselectable=True,
                                           children=[
                                               antenna_card_hover(app,
                                                                  dmc.Chip(s.name, value=s.codename,
                                                                           id={'type': 'antenna-chip',
                                                                               'index': s.codename},
                                                                           color='#004990',
                                                                           persistence=True,
                                                                           styles={'display': 'grid',
                                                                                   'grid-template-columns':
                                                                                   'repeat(auto-fit, '
                                                                                   'minmax(10rem, 1fr))'}),
                                                                  s, show_wavelengths)
                                               for s in observation._STATIONS
                                               if s.codename not in grouped_codenames
                                           ]),
                         ] + grouped_components,
                         className='mb-2 flex',
                         style={'display': 'inline-flex', 'gap': '5px', 'justify-content': 'center',
                                'flex-wrap': 'wrap'})])


def duration() -> html.Div:
    """Returns the DIV with the inputs to set the observation duration (in hours, `duration`)
    and the percentage of time spent on target (`onsourcetime`).

    Returns
    -------
    html.Div
        The duration-and-on-source-time input DIV.
    """
    return html.Div([html.H4("Duration of the Observation", className='text-dark font-weight-bold mb-1'),
                     html.Div(className='col-12', children=[
                        html.Div(className='row d-flex align-items-bottom', children=[
                            html.Div(className='col-5', children=[
                                html.Div(className='row form-group', children=[
                                    html.Label('In hours', htmlFor='duration'),
                                    dbc.Input(id='duration', value=24, type='number', className='form-control',
                                            placeholder="In hours", min=0.1,
                                            persistence=True, debounce=True, inputMode='numeric', max=50, step=0.1),
                                    html.Small(id='error_duration', className='form-text text-muted')])]),
                            html.Div(className='col-7', children=[
                                html.Div(className='row form-group', children=[
                                    html.Label('Percentage of observing time spent on target', htmlFor='onsourcetime'),
                                    dcc.Slider(id='onsourcetime', min=20, max=100, step=10, value=70,
                                                marks={i: str(i) for i in range(20, 101, 10)},
                                                persistence=True, allow_direct_input=False)])])])])])

def source_and_epoch_selection() -> html.Div:
    """Returns the DIV with the epoch selection (start date/time, `startdate`/`starttime`, only
    applied to the computation when `switch-specify-epoch` is on) and the button that opens the
    target-sources modal (`target_sources_modal`), plus the chip display of currently added
    target sources (`target-chips-display`).

    Returns
    -------
    html.Div
        The source-and-epoch selection DIV.
    """
    return html.Div([html.H4("Source  &  Epoch", className='text-dark font-weight-bold mb-1'),
                     html.Div(className='col-12', children=[
                        html.Div(className='row d-flex align-items-bottom', children=[
                            html.Div(className='col-6', children=[
                                html.Div(className='row form-group', children=[
                                    dbc.Switch(label='Specify an epoch', value=False,
                                                id='switch-specify-epoch', persistence=True),
                                    html.Div(id='epoch-selection-div', className='', children=[
                                        html.Div(className='row', children=[
                                            html.Label('Start of observation (UTC)', htmlFor='starttime'),
                                            html.Div(className='row mx-0', children=[
                                                html.Div(className='col-12 px-0 mx-0 d-flex align-items-center gap-2', children=[
                                                    dcc.DatePickerSingle(id='startdate', min_date_allowed=dt(1950, 1, 1),
                                                                         max_date_allowed=dt(2100, 1, 1),
                                                                         initial_visible_month=dt.today(),
                                                                         date=dt.today().date(),
                                                                         className='date-picker',
                                                                         display_format='DD MMM Y'),
                                                    # dcc.Input(id='startdate', type='date',
                                                    #           value=dt.today().date().isoformat(),
                                                    #           min="1950-01-01", max="2100-01-01",
                                                    #           className='date-picker',
                                                    #           display_format='DD MMM Y',
                                                    #           persistence=True, debounce=True),
                                                    dcc.Dropdown(id='starttime', placeholder="Start time (UTC)",
                                                                 value="00:00", clearable=False, searchable=False,
                                                                 options=[{'label': f"{hm//60:02n}:{hm % 60:02n}",
                                                                           'value': f"{hm//60:02n}:{hm % 60:02n}"}
                                                                          for hm in range(0, 24*60, 15)],  # type: ignore[arg-type]
                                                                 persistence=True, className='form-hour')])]),
                                            html.Div(className='row', children=[
                                                html.Small(id='error_starttime', style={'color': 'red'},
                                                        className='form-text text-muted')])])])])]),
                            html.Div(className='col-6', children=[
                                html.Div(className='row form-group', children=[
                                    # html.Label('Target sources', className='form-label fw-bold'),
                                    dbc.Button("Add target sources",
                                               id='button-open-source-modal',
                                               color='primary', outline=True,
                                               className='btn btn-sm w-100'),
                                    html.Div(id='target-chips-display', className='mb-2',
                                             children=html.Small("No target sources added yet.",
                                                                 className='text-muted'))])])])]),
                    ])


def target_sources_modal() -> html.Div:
    """Returns a modal window allowing users to add/remove multiple target sources.

    Sources can be added by typing a name or coordinates, or by uploading a text file
    where each line contains one source. The list of selected target sources is
    persisted in the `store-targets` `dcc.Store`.

    Returns
    -------
    html.Div
        The target-sources modal window.
    """
    return dbc.Modal([
        dbc.ModalHeader(dbc.ModalTitle("Target Sources")),
        dbc.ModalBody([
            html.Div(className='mb-3', children=[
                html.Label("Add a source by name or coordinates:",
                           htmlFor='modal-source-input', className='form-label'),
                dbc.InputGroup([
                    dbc.Input(id='modal-source-input', type='text',
                              placeholder="3C273  or  hh:mm:ss dd:mm:ss",
                              debounce=False, n_submit=0),
                    dbc.Button("Add", id='button-add-source',
                               color='success', n_clicks=0,
                               style={'font-weight': 'bold',
                                      'border-top-left-radius': '0.375rem',
                                      'border-bottom-left-radius': '0.375rem'})]),
                html.Small(id='modal-source-feedback',
                           className='form-text text-muted')]),
            html.Hr(),
            html.Div(className='mb-3', children=[
                html.Label("Or upload a text file (one source per line):",
                           className='form-label'),
                dcc.Upload(id='upload-sources',
                           children=html.Div([
                               html.I(className='fa fa-upload me-2'),
                               'Drag & drop or ',
                               html.A('select a file',
                                      style={'color': '#004990', 'cursor': 'pointer',
                                             'text-decoration': 'underline'})]),
                           className='upload-zone p-3 text-center',
                           style={'border': '2px dashed #004990',
                                  'border-radius': '8px',
                                  'background-color': 'rgba(0, 73, 144, 0.05)'},
                           multiple=False, accept='.txt,.csv,.cat,.lis,.list'),
                html.Small(id='upload-sources-feedback',
                           className='form-text text-muted')]),
            html.Hr(),
            html.Div([
                html.H6("Current target sources", id='modal-source-list-header',
                        className='mb-2'),
                html.Div(id='modal-source-list',
                         children=html.P("No sources added yet.",
                                         className='text-muted text-center my-3'))])]),
        dbc.ModalFooter([
            dbc.Button("Clear all", id='button-clear-sources',
                       color='danger', outline=True, n_clicks=0),
            dbc.Button("Done", id='button-close-source-modal',
                       color='secondary', className='ms-auto', n_clicks=0)])
    ], id='source-modal', is_open=False, size='lg', scrollable=True, backdrop=True)

def correlations() -> html.Div:
    """Creates the inputs to specify if the user wants continuum correlation and/or spectral line.

    Returns
    -------
    html.Div
        The observation-setup input DIV.
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
                                             className='form-dropdown',
                                             options=tuple({'label': drl, 'value': dr}
                                                           for dr, drl in fs.data_rates.items()),  # type: ignore[arg-type]
                                             value=2048, persistence=True, clearable=False, searchable=False),
                                html.Label(id='bandwidth-label', style={'color': '#999999'}, children='',
                                           htmlFor="datarate")]),
                            html.Div(className='col-6', children=[
                                html.Div(className='form-group', children=[
                                    dcc.Dropdown(id='subbands', placeholder="Select no. subbands...",
                                                 className='form-dropdown',
                                                 options=tuple({'label': sbl, 'value': sb}  # type: ignore
                                                               for sb, sbl in fs.subbands.items())[::-1],
                                                 value=8, persistence=True, clearable=False, searchable=False),
                                    dcc.Dropdown(id='channels',
                                                 placeholder="Select no. channels per subband...",
                                                 className='form-dropdown',
                                                 options=tuple({'label': chl, 'value': ch}  # type: ignore
                                                               for ch, chl in fs.channels.items())[::-1],
                                                 value=64, persistence=True, clearable=False, searchable=False),
                                    dcc.Dropdown(id='pols', placeholder="Polarizations...",
                                                 className='form-dropdown',
                                                 options=tuple({'label': pl, 'value': p}  # type: ignore
                                                               for p, pl in fs.polarizations.items()),
                                                 value=4, persistence=True, clearable=False, searchable=False),
                                    dcc.Dropdown(id='inttime', placeholder="Integration time...",
                                                 className='form-dropdown',
                                                 options=tuple({'label': itl, 'value': it}  # type: ignore
                                                               for it, itl in fs.inttimes.items()),
                                                 value=2, persistence=True, clearable=False, searchable=False)])])])])])
