import numpy as np
from astropy import units as u
import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import plotly.express as px


def tooltip(message, idname, trigger='?', placement='right',trigger_is_sup=True,  **kwargs):
    """Defines a tooltip (popover) that will be shown in the rendered page.
    It will place a <sup>`trigger`</sup> in the page within a span tag.
    Returns the list with all html elements.
    """
    if trigger_is_sup:
        return [html.Span(children=html.Sup(trigger, className='popover-link'), id=idname),
            dbc.Tooltip(message, target=idname, placement=placement, **kwargs)]
    else:
        return [html.Span(children=trigger, className='popover-link', id=idname),
            dbc.Tooltip(message, target=idname, placement=placement, **kwargs)]


def tooltip_card(a_card, idname, trigger, placement='right', **kwargs):
    """Defines a tooltip (popover) that shows a card.
    """
    return [html.Span(children=trigger, className='popover-link', id=idname),
            dbc.Tooltip(a_card, target=idname, placement=placement,
                innerClassName='tooltip-card-inner', hide_arrow=True,
                **kwargs)]


def create_accordion_card(title, text, id, is_open=True):
    """Given a title (header) and a text (which can be either text, a dcc/html object),
    it will return a dbc.Card object which is one input for an accordion.
    """
    card_header = dbc.CardHeader(html.H2(dbc.Button(title, color='link',
                    id=f"group-{id}-toggle", className='')), className='accordion-header')
    card_body = dbc.Collapse(dbc.CardBody(text), id=f"collapse-{id}", is_open=is_open,
                    className='accordion-collapse')
    return dbc.Card([card_header, card_body], className='accordion-card')


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




#################################################################################
# It is all about cards in this section

def antenna_card(app, station):
    """Generates a card showing the basic information for the given station
    """
    s = lambda st : st[::-1].replace(' ,',' dna ',1)[::-1]
    card = dbc.Card([
        dbc.CardImg(src=app.get_asset_url(f"ant-{station.name.replace(' ','_').lower()}.jpg"),
                    top=True, className='card-img'),
        dbc.CardBody([
            html.H4(station.name, className='card-title'),
            html.H6(station.fullname if station.fullname != station.name else '',
                    className='card-title2'),
            html.H6([html.Span(station.country),
                     html.Span(station.diameter, style={'float': 'right'}),
                    ], className='card-subtitle'),
            # html.P(f"&#127462; Participates in {station.all_networks}.\n"
            dcc.Markdown(f"Listed for the {s(station.all_networks)}.\n" if \
                          station.all_networks != '' else '', className='card-text'),
            dcc.Markdown("Can observe at "
                         f"{', '.join([i.replace('cm', '') for i in station.bands])} cm.",
                         className='card-text')
            ])
        ], className='card-antenna')
    return card


def antenna_cards(app, stations):
    cards = dbc.Row([antenna_card(app, s) for s in stations],
            className='row justify-content-center')
    return cards



def summary_card_antennas(app, obs):
    """Generates the summary card with the information about which
    antennas can observe the given source, and the longest/shortest baselines.
    """
    ants_up = obs.is_visible()
    ant_no_obs = []
    for an_ant in ants_up:
        if len(ants_up[an_ant][0]) == 0:
            ant_no_obs.append(an_ant)

    ant_text = ', '.join([ant for ant in obs.stations.codenames if ant not in ant_no_obs]) + '.'
    ant_text = []
    for ant in obs.stations.codenames:
        if ant not in ant_no_obs:
            ant_text += tooltip_card(antenna_card(app, obs.stations[ant]),
                            idname=f"basel-{ant}", trigger=ant, placement='top')
            ant_text += [html.Span(", ")]

    # Remove the trailing ,
    if len(ant_text) > 0:
        ant_text[-1] = html.Span(".")
    else:
        ant_text = ["None."]
    # TODO: This is the worldmap... I think it does not fit here.
    # temp_msg = [ge.worldmap_plot([obs.stations[a] for a in obs.stations.keys() \
    #             if a not in ant_no_obs])]
    temp_msg = []
    temp_msg += [[f"{len(ants_up)-len(ant_no_obs)} participating antennas: ", *ant_text]]
    if len(ant_no_obs) > 0:
        ant_text = []
        for ant in ant_no_obs:
            ant_text += tooltip_card(antenna_card(app, obs.stations[ant]),
                            idname=f"basel-{ant}", trigger=ant, placement='top')
            ant_text += [html.Span(", ")]

        temp_msg += [html.P(className='text-danger', children=["Note that ", *ant_text[:-1],
                                                    " cannot observe the source."])]

    longest_bl = obs.longest_baseline()
    ant_l1, ant_l2 = longest_bl[0].split('-')
    # Using dummy units to allow the conversion
    longest_bl_lambda = longest_bl[1]/obs.wavelength
    longest_bl_lambda = optimal_units(longest_bl_lambda*u.m, [u.Gm, u.Mm, u.km])
    temp_msg += [[*baseline_img(app, is_long=True),
                *tooltip_card(antenna_card(app, obs.stations[ant_l1]), idname='basel-l1',
                             trigger=ant_l1, placement='top'),
                html.Span("-"),
                *tooltip_card(antenna_card(app, obs.stations[ant_l2]), idname='basel-l2',
                             trigger=ant_l2, placement='top'),
                f" is the longest (projected) baseline with {optimal_units(longest_bl[1], [u.km, u.m]):.5n} ({longest_bl_lambda.value:.3n} {longest_bl_lambda.unit.name[0]}\u03BB)."]]

    shortest_bl = obs.shortest_baseline()
    ant_s1, ant_s2 = shortest_bl[0].split('-')
    # Using dummy units to allow the conversion
    shortest_bl_lambda = shortest_bl[1]/obs.wavelength
    shortest_bl_lambda = optimal_units(shortest_bl_lambda*u.m, [u.Gm, u.Mm, u.km])
    temp_msg += [[*baseline_img(app, is_long=False),
                *tooltip_card(antenna_card(app, obs.stations[ant_s1]), idname='basel-s1',
                             trigger=ant_s1, placement='top'),
                html.Span("-"),
                *tooltip_card(antenna_card(app, obs.stations[ant_s2]), idname='basel-s2',
                             trigger=ant_s2, placement='top'),
                 f" is the shortest one with {optimal_units(shortest_bl[1], [u.km, u.m]):.5n} ({shortest_bl_lambda.value:.3n} {shortest_bl_lambda.unit.name[0]}\u03BB)."]]
    # )]
    return create_sensitivity_card('Antennas', temp_msg)


def summary_card_beam(app, obs):
    """Creates a summary card showing the expected synthesized beam.
    """
    synthbeam = obs.synthesized_beam()
    synthbeam_units = optimal_units(synthbeam['bmaj'], [u.arcsec, u.mas, u.uas]).unit
    temp_msg = [html.Div(className='row', style={'height': '1rem'}),
                html.Div(className='row justify-content-center',
              children=[ellipse(bmaj="5rem",
                        bmin=f"{5*synthbeam['bmin'].to(u.mas)/synthbeam['bmaj'].to(u.mas)}rem",
                        pa=f"{synthbeam['pa'].to(u.deg).value}deg")])]
    # TODO: Check that the rotation is the correct.
    temp_msg += [html.Br(), html.P([f"The expected synthesized beam will be approx. {synthbeam['bmaj'].to(synthbeam_units).value:.3n} x {synthbeam['bmin'].to(synthbeam_units):.3n}", html.Sup("2"), \
            f", PA = {synthbeam['pa']:.3n}."])]
    temp_msg += [html.P("Note that the synthesized beam can significantly change depending "
                        "on the weighting used during imaging.")]
    return create_sensitivity_card('Resolution', temp_msg)

def summary_card_times(app, obs):
    """Creates a summary card showing the observing times, and the resulting data size.
    """
    prtobstimes = obs.print_obs_times()
    if '\n' in prtobstimes:
        tmp = [html.Span(t) for t in obs.print_obs_times().split('\n')]
        for i in range(len(tmp)-1):
            tmp.insert(2*i+1, html.Br())
            temp_msg = [tmp]
    else:
        temp_msg = [f"{obs.print_obs_times()}."]

    temp_msg += [f"The observation lasts for {optimal_units(obs.duration, [u.h, u.min, u.s, u.ms]):.3n}, of which {optimal_units(obs.ontarget_time, [u.h, u.min, u.s, u.ms]):.3n} are on target."]
    n_files = int(np.ceil(obs.datasize()/(2.0*u.GB)))
    if n_files < 10:
        img_fits = [html.Img(src=app.get_asset_url("icon-file.svg"), height='35rem',
                    style={'display': 'inline'}) for i in range(n_files)]
    else:
        img_fits = [html.Img(src=app.get_asset_url("icon-file.svg"), height='35rem',
                    style={'display': 'inline'}) for i in range(10)]
        img_fits += ["+"]

    temp_msg += [[*img_fits, html.Br(),
        f"With a time integration of {optimal_units(obs.inttime, [u.s,u.ms,u.us]):.2n} the "
        f"expected FITS file size is "
        f"{optimal_units(obs.datasize(), [u.TB, u.GB, u.MB, u.kB]):.3n} "
        f"(divided in {n_files} 2-GB files)."]]
    return create_sensitivity_card('Observing Time', temp_msg)


def summary_card_frequency(app, obs):
    pol_dict = {1: 'single', 2: 'dual', 4: 'full'}
    bw = optimal_units(obs.bandwidth, [u.GHz, u.MHz, u.kHz])
    bwwl = optimal_units((30*u.cm/(obs.bandwidth.to(u.GHz).value)), [u.m, u.cm, u.mm])
    temp_msg = [html.Div(className='row justify-content-center',
                children=[html.Div(className='bandpass', style={'height': '4rem',
                'width': f"{90/obs.subbands}%"})]*obs.subbands)]
    temp_msg += [html.Div(className='row justify-content-center',
                 children=[html.Table(className='baseline-table',
                           style={'width': '100%', 'font-size': '0.8rem', 'margin-top':'-1rem'},
                           children=[html.Tr([
                                html.Td(f"-{bw/2:.3n}", style={'text-align': 'left'}),
                                html.Td(f"{optimal_units(obs.frequency, [u.GHz, u.MHz]):.4n}",
                                        style={'text-align': 'center'}),
                                html.Td(f"+{bw/2:.3n}", style={'text-align': 'right'})
                            ])]
                )])]
    temp_msg += [f"The central frequency is {optimal_units(obs.frequency, [u.GHz, u.MHz]):.3n} ({optimal_units(obs.wavelength, [u.m, u.cm, u.mm]):.2n})."]
    temp_msg += [f"The total bandwidth of {optimal_units(obs.bandwidth, [u.GHz, u.MHz, u.kHz]):.3n} will be divided into {obs.subbands} subbands of {optimal_units(obs.bandwidth/obs.subbands, [u.GHz, u.MHz, u.kHz]):.3n} each, with {obs.channels} channels ({optimal_units(obs.bandwidth/(obs.subbands*obs.channels), [u.GHz, u.MHz, u.kHz, u.Hz]):.3n} wide)."]
    temp_msg += [f"Recording {pol_dict[obs.polarizations]} circular polarization."]
    return create_sensitivity_card('Frequency Setup', temp_msg)



def summary_card_fov(app, obs):
    """Creates a summary card showing the expected FoV limitations.
    """
    # primary_beam =
    shortest_bl = obs.shortest_baseline()[1]
    largest_ang_scales = ((2.063e8*u.mas)*(obs.wavelength/shortest_bl)).to(u.mas)
    pb_scale = ((2.063e8*u.mas)*(obs.wavelength/(100*u.m))).to(u.arcsec)
    bw_smearing = obs.bandwidth_smearing()
    tm_smearing = obs.time_smearing()
    smearing_ratio = bw_smearing/tm_smearing
    smearing_ratio = smearing_ratio if smearing_ratio <= 1.0 else 1/smearing_ratio

    temp_msg = [html.Div(className='row', style={'height': '1rem'}),
            html.Div(className='row justify-content-center',
                    children=[ellipse(bmaj="5rem", bmin="5rem", pa="0deg"),
                      ellipse(bmaj=f"2rem", bmin=f"2rem", pa="0deg", color="#F0959B",
                              z_index=3, position='absolute', margin_top='7%'),
                      ellipse(bmaj=f"{2*smearing_ratio}rem",
                              bmin=f"{2*smearing_ratio}rem",
                              pa="0deg", color="white", z_index=4, position='absolute',
                              margin_top=f'{8+2*bw_smearing/tm_smearing}%')])]
    temp_msg += [f"The Field of View would be limited by time smearing to {optimal_units(tm_smearing, [u.arcmin, u.arcsec]):.3n} and by frequency smearing to {optimal_units(bw_smearing, [u.arcmin, u.arcsec]):.3n} (considering a 10% loss)."]
    temp_msg += [f"Considering the shortest baseline in the array, "
    "you will filter out emission on angular scales larger than "
    f"{optimal_units(largest_ang_scales, [u.arcmin, u.arcsec, u.mas]):.3n}."]
    return create_sensitivity_card('FoV limitations', temp_msg)



def summary_card_rms(app, obs):
    """Creates a summary card showing the reached sensitivity.
    """
    rms = optimal_units(obs.thermal_noise(), [u.MJy, u.kJy, u.Jy, u.mJy, u.uJy])
    rms_time = optimal_units(obs.thermal_noise()/np.sqrt(obs.inttime/obs.duration),
                             [u.MJy, u.kJy, u.Jy, u.mJy, u.uJy])
    rms_channel = optimal_units(rms*np.sqrt(obs.subbands*obs.channels),
                                [u.MJy, u.kJy, u.Jy, u.mJy, u.uJy])
    temp_msg = [html.Div(className='row', style={'height': '0.7rem'}),
                html.Div(className='row justify-content-center',
                children=html.Img(src=app.get_asset_url("waves.svg"), width='100%',
                                  height='75rem', style={'display': 'inline'}))]
    temp_msg += [html.P(f"Considering the sensitivities of the antennas, "
                        f"the estimated thermal noise is {rms:.3n}/beam.")]
    temp_msg += [html.P(f"This would imply a rms of {rms_channel:.3n}/beam per spectral "
                        f"channel, or approx. {rms_time:.3n}/beam per time integration "
                        f"({optimal_units(obs.inttime, [u.s,u.ms,u.us]):.3n}).")]
    return create_sensitivity_card('Sensitivity', temp_msg)



#################################################################################
# Some small graphical elements


def ellipse(bmaj, bmin, pa, color='#a01d26', z_index=1, position='relative', margin_top=''):
    """Returns a html.Div element that draws an ellipse with a semimajor axis bmaj,
    semiminor axis bmin, and position angle (as defined in radio astronomy) pa.
    bmaj,bmin, pa must be strings recognized by HTML/CSS.
    """
    return html.Div(children=[], style={'width': bmaj, 'height': bmin,
                    'border-radius': '50%', 'background': color, 'position': position,
                    'transform': f"rotate({pa})", 'z-index': z_index,
                    'vertical-align': 'middle', 'margin-top': margin_top})


def baseline_img(app, is_long=True):
    """Returns a HTML element visually representing a baseline.
    The only parameter is a bool that defines if the representation
    is for a long baseline (True) or for a short baseline (False).
    """
    if is_long:
        baseline = html.Td(className='baseline-td-hr',children=html.Hr(className='hr-baseline'),
                            style={'width': '80%', 'padding': '0', 'margin': '0'})
    else:
        baseline = html.Td(className='baseline-td-hr',children=html.Hr(className='hr-baseline'),
                            style={'width': '30%', 'padding': '0', 'margin': '0'})

    return [html.Table(className='baseline-table', children=[
            html.Tr(className='baseline-tr', children=[
                html.Td(className='baseline-td-img', children=html.Img(width='17rem',
                        src=app.get_asset_url("icon-32.png"), alt='Antenna logo',
                        className='d-inline-block align-right img-baseline')),
                baseline,
                html.Td(className='baseline-td-img', children=html.Img(width='17rem',
                        src=app.get_asset_url("icon-32.png"), alt='Antenna logo',
                        className='d-inline-block align-left img-baseline')),
            ])])]




#################################################################################
# Some small graphical elements

def worldmap_plot(antennas):
    data = {"lat": [], "lon": [], "color": [], "symbol": [], "mode": "markers",
            "name": [], "hovertext": []}
    for ant in antennas:
        data["lat"].append(ant.location.lat.value)
        data["lon"].append(ant.location.lon.value)
        data["color"].append('green')
        data["symbol"].append(0.1)
        data["name"].append(ant.name)
        data["hovertext"].append(ant.name)

    fig = px.scatter_geo(data,
                    lat="lat", lon="lon", color="color", text="hovertext",
                    hover_name="name")
    fig.update_layout(autosize=True, hovermode='closest', showlegend=False,
                      margin={'l': 0, 't': 0, 'b': 0, 'r': 0})

    return html.Div(className='row justify-content-center', children=dcc.Graph(figure=fig,
            responsive=True, className='fig-on-card',
                      config={'frameMargins': -10, 'showLink': False, 'displaylogo': False}))


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




