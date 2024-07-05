import numpy as np
from datetime import datetime as dt
import enum
from typing import Optional, Union
from astropy import units as u
from fpdf import FPDF
from dash import dcc
from dash import html
import dash_bootstrap_components as dbc
import plotly.express as px
from vlbiplanobs import freqsetups as fs
from vlbiplanobs import stations
from vlbiplanobs import observation




class SourceEpoch(enum.Enum):
    """Enum-type class with the three possible types of observations to set,
    in terms of how to specify the source and epoch to observe:
    - UNSPECIFIED: if neither source nor epoch have been specified.
    - ONLY_SOURCE: only source is specified, not the epoch.
    - SOURCE_AND_EPOCH: both the epoch and the duration are specified.
    """
    UNSPECIFIED = 0
    ONLY_SOURCE = 1
    SOURCE_AND_EPOCH = 2



def tooltip(message: str, idname: str, trigger: str = '?', placement: str = 'right',
            trigger_is_sup: bool = True,  **kwargs) -> list:
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


def tooltip_card(a_card: dbc.Card, idname: str, trigger: str, placement: str = 'right', **kwargs) -> list:
    """Defines a tooltip (popover) that shows a card.
    """
    return [html.Span(children=trigger, className='popover-link', id=idname),
            dbc.Tooltip(a_card, target=idname, placement=placement,
                innerClassName='tooltip-card-inner', hide_arrow=True,
                **kwargs)]


def create_accordion_card(title: str, text: str, id: str, is_open: bool = True) -> dbc.Card:
    """Given a title (header) and a text (which can be either text, a dcc/html object),
    it will return a dbc.Card object which is one input for an accordion.
    """
    card_header = dbc.CardHeader(html.H2(dbc.Button(title, color='link',
                    id=f"group-{id}-toggle", className='')), className='accordion-header')
    card_body = dbc.Collapse(dbc.CardBody(text), id=f"collapse-{id}", is_open=is_open,
                    className='accordion-collapse')
    return dbc.Card([card_header, card_body], className='accordion-card')


def create_sensitivity_card(title: str, message: Union[str, list]) -> list:
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

    # return [html.Div(dbc.Card(
    return [dbc.Card(
            dbc.CardBody([
                html.H5(className='card-title', children=title),
                *ps
                ]), className='card-summary shadow-1-strong')
            ]
    # return [html.Div(className='card-summary m-3', children=[
    #         html.Div(className='card-body', children=[
    #             html.H5(className='card-title', children=title)] + ps)])
    #     ]




#################################################################################
# It is all about cards in this section

def antenna_card(app, station: stations.Station) -> dbc.Card:
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
            dcc.Markdown([f"Listed for the {s(station.all_networks)}.\n" if \
                          station.all_networks != '' else '', "Can observe at "
                         f"{', '.join([i.replace('cm', '') for i in station.bands])} cm."],
                         className='card-text')
            ])
        ], className='card-antenna')
    return card



def antenna_cards(app, stations: list) -> dbc.Row:
    cards = dbc.Row([antenna_card(app, s) for s in stations],
            className='row justify-content-center')
    return cards



def network_card(app, network_acr: str, network_name: str, body,
                 network_img: Optional[str] = None) -> dbc.Card:
    """Generates a card showing the basic information from a network.
    """
    card = dbc.Card([
        dbc.CardImg(src=app.get_asset_url(network_img if network_img is not None else \
                                          f"network-{network_acr.lower().replace(' ', '_')}.png"),
                    top=True, className='card-img'),
        dbc.CardBody([
            html.H4(network_acr, className='card-title'),
            html.H6(network_name, className='card-title2'),
            *[dbc.Checklist(id='initial-e-EVN', className='checkbox', persistence=True,
                    options=[{'label': 'e-EVN (real time) mode', 'value': False}]) \
              if network_acr == 'EVN' else html.Span()],
            ]),
        dbc.CardFooter(
            dbc.Checklist(id=f'network-{network_acr.lower()}', className='checkbox', persistence=True,
                          options=[{'label': ' Select network',
                                    'value': True}], value=[]))
        ], className="card-network col-sm-3 m-4 shadow-1-strong")
    return card



def summary_card_worldmap(app, obs: observation.Observation) -> html.Div:
    """Generates a world map with the participating stations
    """
    ants_up = obs.is_visible()
    ant_no_obs = []
    for an_ant in ants_up:
        if len(ants_up[an_ant][0]) == 0:
            ant_no_obs.append(an_ant)

    return worldmap_plot([obs.stations[a] for a in obs.stations.codenames \
                if a not in ant_no_obs])



def summary_card_antennas(app, obs: observation.Observation) -> list:
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
    temp_msg = []
    temp_msg += [[f"{len(ants_up)-len(ant_no_obs)} participating antennas: ", *ant_text]]
    if len(ant_no_obs) > 0:
        ant_text = []
        for ant in ant_no_obs:
            ant_text += tooltip_card(antenna_card(app, obs.stations[ant]),
                            idname=f"basel-{ant}", trigger=ant, placement='top')
            ant_text += [html.Span(", ")]

        temp_msg += [html.P(className='text-danger', children=["Note that ", *ant_text[:-1],
                            " cannot observe the source during the planned observation."])]

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
                        pa=f"{-synthbeam['pa'].to(u.deg).value+90}deg")])]
    temp_msg += [html.Br(), html.P([f"The expected synthesized beam will be approx. {synthbeam['bmaj'].to(synthbeam_units).value:.3n} x {synthbeam['bmin'].to(synthbeam_units):.3n}", html.Sup("2"), \
            f", PA = {synthbeam['pa']:.3n}."])]
    temp_msg += [html.P("Note that the synthesized beam can significantly change depending "
                        "on the weighting used during imaging (natural weighting assumed here).")]
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
        f"{optimal_units(obs.datasize(), [u.TB, u.GB, u.MB, u.kB]):.3n}. "]]
        # f"(divided in {n_files} 2-GB files)."]]
    return create_sensitivity_card('Observing Time', temp_msg)


def summary_card_frequency(app, obs):
    pol_dict = {1: 'single', 2: 'dual', 4: 'full'}
    bw = optimal_units(obs.bandwidth, [u.GHz, u.MHz, u.kHz])
    vel = optimal_units((2.9979e5*u.km/u.s)*(obs.bandwidth/(2*obs.frequency)).decompose(), [u.km/u.s, u.m/u.s])
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
                            ]),
                            html.Tr([
                                html.Td(f"+{vel:.5n}", style={'text-align': 'left'}),
                                html.Td("0", style={'text-align': 'center'}),
                                html.Td(f"-{vel:.5n}", style={'text-align': 'right'})
                            ])]
                )])]
    temp_msg += [f"The central frequency is {optimal_units(obs.frequency, [u.GHz, u.MHz]):.3n} ({optimal_units(obs.wavelength, [u.m, u.cm, u.mm]):.2n})."]
    if obs.subbands > 1:
        temp_msg += [f"The total bandwidth of {optimal_units(obs.bandwidth, [u.GHz, u.MHz, u.kHz]):.3n} will be divided into {obs.subbands} subbands of {optimal_units(obs.bandwidth/obs.subbands, [u.GHz, u.MHz, u.kHz]):.3n} each, with {obs.channels} channels ({optimal_units(obs.bandwidth/(obs.subbands*obs.channels), [u.GHz, u.MHz, u.kHz, u.Hz]):.3n}, or {optimal_units(2*vel/(obs.subbands*obs.channels), [u.km/u.s, u.m/u.s]):.3n}, wide)."]
    else:
        temp_msg += [f"The total bandwidth of {optimal_units(obs.bandwidth, [u.GHz, u.MHz, u.kHz]):.3n} will be recorded in {obs.subbands} subband, with {obs.channels} channels ({optimal_units(obs.bandwidth/(obs.subbands*obs.channels), [u.GHz, u.MHz, u.kHz, u.Hz]):.3n}, or {optimal_units(2*vel/(obs.subbands*obs.channels), [u.km/u.s, u.m/u.s]):.3n}, wide)."]
    temp_msg += [f"Recording {pol_dict[obs.polarizations]} circular polarization."]
    return create_sensitivity_card('Frequency Setup', temp_msg)



def summary_card_fov(app, obs):
    """Creates a summary card showing the expected FoV limitations.
    """
    shortest_bl = obs.shortest_baseline()[1]
    largest_ang_scales = ((2.063e8*u.mas)*(obs.wavelength/shortest_bl)).to(u.mas)
    pb_scale = ((2.063e8*u.mas)*(obs.wavelength/(100*u.m))).to(u.arcsec)
    bw_smearing = obs.bandwidth_smearing()
    tm_smearing = obs.time_smearing()
    smearing_ratio = bw_smearing/tm_smearing
    smearing_ratio = smearing_ratio if smearing_ratio <= 1.0 else 1/smearing_ratio
    temp_msg = [html.Div(className='row', style={'height': '1rem'}),
            html.Div(className='row justify-content-center align-self-center',
                    children=[ellipse(bmaj="5rem", bmin="5rem", pa="0deg"),
                      ellipse(bmaj=f"2rem", bmin=f"2rem", pa="0deg", color="#F0959B",
                              z_index=3, position='absolute', className='align-self-center'),
                      ellipse(bmaj=f"{2*smearing_ratio}rem",
                              bmin=f"{2*smearing_ratio}rem",
                              pa="0deg", color="white", z_index=4, position='absolute',
                              className='align-self-center')])]
    temp_msg += [f"The Field of View would be limited by time smearing to a radius of "
                 f"{optimal_units(tm_smearing, [u.arcmin, u.arcsec]):.3n} and by frequency smearing to "
                 f"{optimal_units(bw_smearing, [u.arcmin, u.arcsec]):.3n} (considering a 10% loss), "
                 "if no further time/frequency averaging is performed."]
    temp_msg += [f"Considering the shortest baseline in the array, "
                 "you will resolve out emission on angular scales larger than "
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
                children=html.Img(src=app.get_asset_url("waves.png"), width='100%',
                                  height='75rem', style={'display': 'inline'}))]
    temp_msg += [html.P(f"The expected rms thermal noise for your target is {rms:.3n}/beam using natural weighting during imaging. Note that ~20% higher values may be expected for RFI-contaminated bands.")]
    temp_msg += [html.P(f"The achieved sensitivity implies a rms of {rms_channel:.3n}/beam per spectral "
                        f"channel, or approx. {rms_time:.3n}/beam per time integration "
                        f"({optimal_units(obs.inttime, [u.s,u.ms,u.us]):.3n}).")]
    return create_sensitivity_card('Sensitivity', temp_msg)



def summary_printable(app, obs):
    """Returns a html.Div object with a HTML format.
    """
    rms = optimal_units(obs.thermal_noise(), [u.MJy, u.kJy, u.Jy, u.mJy, u.uJy])
    rms_channel = optimal_units(rms*np.sqrt(obs.subbands*obs.channels),
                        [u.MJy, u.kJy, u.Jy, u.mJy, u.uJy])
    bw_smearing = obs.bandwidth_smearing()
    tm_smearing = obs.time_smearing()
    pol_dict = {1: 'single', 2: 'dual', 4: 'full'}
    synthbeam = obs.synthesized_beam()
    synthbeam_units = optimal_units(synthbeam['bmaj'], [u.arcsec, u.mas, u.uas]).unit
    ants_up = obs.is_visible()
    ant_no_obs = []
    for an_ant in ants_up:
        if len(ants_up[an_ant][0]) == 0:
            ant_no_obs.append(an_ant)

    ants = [ant for ant in obs.stations.codenames if ant not in ant_no_obs]
    if obs.target is not None:
        target = f"Target source: {obs.target.coord.to_string('hmsdms')}{' ('+obs.target.name+').' if obs.target.name != 'Source' else '.'}"
    else:
        target = "Target source: unspecified."

    epoch = f"{obs.print_obs_times()}" if obs._fixed_time else \
        f"Observing epoch: unspecified.\nTotal observing time: {optimal_units(obs.duration, [u.h, u.min, u.s]):.3n}"
    # Create the PDF and returns it
    pdf = FPDF(orientation='P', unit='mm', format='A4')
    pdf.compress = False
    pdf_w = 210
    pdf_h = 297
    # self.image(sctplt,  link='', type='', w=1586/80, h=1920/80)
    pdf.add_page()
    pdf.ln(1)
    pdf.set_xy(0.0, 0.0)
    pdf.set_font('Arial', 'B', 16)
    # self.set_text_color(220, 50, 50)
    pdf.cell(w=210.0, h=40.0, align='C', txt="EVN Observation Planner - Summary Report", border=0)
    pdf.set_xy(10.0, 40.0)
    pdf.set_font('Arial', 'B', 12)
    pdf.multi_cell(0, 0, "Schedule")
    pdf.set_xy(10.0, 45.0)
    pdf.set_font('Arial', '', 10)
    pdf.multi_cell(0, 6, f"{epoch}.\n" \
       f"{optimal_units(obs.ontarget_time, [u.h, u.min, u.s, u.ms]):.3n} are on target.\n" \
       f"{target}\n" \
       f"Output FITS file size: {optimal_units(obs.datasize(), [u.TB, u.GB, u.MB, u.kB]):.3n}.\n\n\n")
    pdf.set_xy(10.0, 85.0)
    pdf.set_font('Arial', 'B', 12)
    pdf.multi_cell(0, 0, "Frequency Setup")
    pdf.set_xy(10.0, 90.0)
    pdf.set_font('Arial', '', 10)
    pdf.multi_cell(0, 6, f"Central frequency: {optimal_units(obs.frequency, [u.GHz, u.MHz]):.3n} " \
       f"({optimal_units(obs.wavelength, [u.m, u.cm, u.mm]):.2n}).\n" \
       f"{obs.subbands} subbands of {optimal_units(obs.bandwidth/obs.subbands, [u.GHz, u.MHz, u.kHz]):.3n} each.\n" \
       f"Channels per subband: {obs.channels}.\n" \
       f"Polarization: {pol_dict[obs.polarizations]}.\n" \
       f"time integration: {optimal_units(obs.inttime, [u.s,u.ms,u.us]):.2n}.\n\n\n")
    pdf.set_xy(10.0, 135.0)
    pdf.set_font('Arial', 'B', 12)
    pdf.multi_cell(0, 0, "VLBI Network")
    pdf.set_xy(10.0, 140.0)
    pdf.set_font('Arial', '', 10)
    pdf.multi_cell(0, 6, f"Participating antennas: {', '.join(ants)}.\n" \
       f"The expected synthesized beam will be approx. {synthbeam['bmaj'].to(synthbeam_units).value:.3n} x {synthbeam['bmin'].to(synthbeam_units):.3n}, PA = {synthbeam['pa']:.3n}.\n" \
       f"Expected rms thermal noise level: {rms:.3n}/beam.\n" \
       f"Per spectral channel: {rms_channel:.3n}/beam.\n" \
       f"Time smearing (10% loss): {optimal_units(tm_smearing, [u.arcmin, u.arcsec]):.3n}.\n" \
       f"Frequency smearing (10% loss): {optimal_units(bw_smearing, [u.arcmin, u.arcsec]):.3n}.\n")

    return pdf
    # Single string version
    return "EVN Observation Planner - Summary Report\n\n\n" \
           "- Schedule\n\n" \
           f"{obs.print_obs_times()}.\n" \
           f"{optimal_units(obs.ontarget_time, [u.h, u.min, u.s, u.ms]):.3n} are on target.\n" \
           f"{target}\n" \
           f"FITS file size: {optimal_units(obs.datasize(), [u.TB, u.GB, u.MB, u.kB]):.3n}.\n\n\n" \
           "- Frequency Setup\n\n" \
           f"Central frequency: {optimal_units(obs.frequency, [u.GHz, u.MHz]):.3n} " \
           f"({optimal_units(obs.wavelength, [u.m, u.cm, u.mm]):.2n}).\n" \
           f"{obs.subbands} subbands of {optimal_units(obs.bandwidth/obs.subbands, [u.GHz, u.MHz, u.kHz]):.3n} each.\n" \
           f"Channels per subband: {obs.channels}.\n" \
           f"Polarization: {pol_dict[obs.polarizations]}.\n" \
           f"time integration: {optimal_units(obs.inttime, [u.s,u.ms,u.us]):.2n}.\n\n\n" \
           "- VLBI Network\n\n" \
           f"Participating antennas: {ants}.\n" \
           f"The expected synthesized beam will be approx. {synthbeam['bmaj'].to(synthbeam_units).value:.3n} " \
           f"x {synthbeam['bmin'].to(synthbeam_units):.3n}, PA = {synthbeam['pa']:.3n}.\n" \
           f"Expected rms thermal noise level: {rms:.3n}/beam.\n" \
           f"Per spectral channel: {rms_channel:.3n}/beam.\n" \
           f"Time smearing (10% loss): {optimal_units(tm_smearing, [u.arcmin, u.arcsec]):.3n}\n" \
           f"Frequency smearing (10% loss): {optimal_units(bw_smearing, [u.arcmin, u.arcsec]):.3n}\n" \
           f"Per spectral channel: {rms_channel:.3n}/beam.\n"
    # Dash version
    return html.Div([
            html.H3("EVN Observation Planner - Report"),
            html.H5("Schedule"),
            html.P(f"{obs.print_obs_times()}."),
            html.P(f"{optimal_units(obs.ontarget_time, [u.h, u.min, u.s, u.ms]):.3n} are on target."),
            html.P(f"Target source: {obs.target.coord.to_string('hmsdms')}"
                   f"{' ('+obs.target.name+').' if obs.target.name != 'Source' else '.'}"
                   if obs.target is not None else 'Target source: Unspecified.'),
            html.P(f"FITS file size: {optimal_units(obs.datasize(), [u.TB, u.GB, u.MB, u.kB]):.3n}."),

            html.H5("Frequency Setup"),
            html.P(f"Central frequency: {optimal_units(obs.frequency, [u.GHz, u.MHz]):.3n} "
                   f"({optimal_units(obs.wavelength, [u.m, u.cm, u.mm]):.2n})."),
            html.P(f"{obs.subbands} subbands of "
                   f"{optimal_units(obs.bandwidth/obs.subbands, [u.GHz, u.MHz, u.kHz]):.3n} each."),
            html.P(f"Channels per subband: {obs.channels}."),
            html.P(f"Polarization: {pol_dict[obs.polarizations]}."),
            html.P(f"time integration: {optimal_units(obs.inttime, [u.s,u.ms,u.us]):.2n}."),
            html.P(f""),

            html.H5("VLBI Network"),
            html.P(f"Participating antennas: {ants}."),
            html.P("The expected synthesized beam will be approx. "
                   f"{synthbeam['bmaj'].to(synthbeam_units).value:.3n} x "
                   f"{synthbeam['bmin'].to(synthbeam_units):.3n}, PA = {synthbeam['pa']:.3n}."),
            html.P(f"Expected rms thermal noise level: {rms:.3n}/beam."),
            html.P(f"Per spectral channel: {rms_channel:.3n}/beam."),
            html.P(f"Time smearing (10% loss): {optimal_units(tm_smearing, [u.arcmin, u.arcsec]):.3n}"),
            html.P(f"Frequency smearing (10% loss): {optimal_units(bw_smearing, [u.arcmin, u.arcsec]):.3n}"),
            html.P(f"Per spectral channel: {rms_channel:.3n}/beam."),
            # TODO: Add figure elevations per antenna
        ])


#################################################################################
# Some small graphical elements


def ellipse(bmaj, bmin, pa, color='#a01d26', z_index=1, position='relative', margin_top='', className=''):
    """Returns a html.Div element that draws an ellipse with a semimajor axis bmaj,
    semiminor axis bmin, and position angle (as defined in radio astronomy) pa.
    bmaj,bmin, pa must be strings recognized by HTML/CSS.
    """
    return html.Div(children=[], style={'width': bmaj, 'height': bmin,
                    'border-radius': '50%', 'background': color, 'position': position,
                    'transform': f"rotate({pa})", 'z-index': z_index,
                    'vertical-align': 'middle', 'margin-top': margin_top}, className=className)


def baseline_img(app, is_long=True):
    """Returns a HTML element visually representing a baseline.
    The only parameter is a bool that defines if the representation
    is for a long baseline (True) or for a short baseline (False).
    """
    if is_long:
        baseline = html.Td(className='baseline-td-hr', children=html.Hr(className='hr-baseline'),
                           style={'width': '80%', 'padding': '0', 'margin': '0'})
    else:
        baseline = html.Td(className='baseline-td-hr', children=html.Hr(className='hr-baseline'),
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
            "name": [], "text": [], "hovertemplate": []}
    for ant in antennas:
        data["lat"].append(ant.location.lat.value)
        data["lon"].append(ant.location.lon.value)
        data["color"].append('green')
        data["symbol"].append(0.1)
        data["name"].append(ant.name)
        data["text"].append(f"{ant.name}<br>({ant.country})<br> {ant.diameter}")
        data["hovertemplate"].append(f"{ant.name}<br>({ant.country})<br> {ant.diameter}")

    fig = px.scatter_geo(data,
                    lat="lat", lon="lon", text="name",
                    hover_name="text", hover_data=None)
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


#################################################################################
# Cards from the initial window where user keeps selecting options


def initial_window_start(app):
    """Very first window where the user can choose between going to the main window
    directly to manually selecting everything or going through the "wizard" guided mode.
    """
    return [
            html.Div(className='row card-deck justify-content-center text-center', children=[
                html.Button(id='button-initial-wizard', className='card-button btn btn-gray m-4',
                    title='Helps you to configure your observation in a step-by-step process.',
                    children=dbc.Card(dbc.CardBody([
                        html.H5("Guided mode", className='card-title'),
                        html.Img(height='80rem', src=app.get_asset_url('icon-newbie.png'),
                                 className='card-text m-3'),
                        html.Br(),
                        html.P("Easy step-by-step setup", className='card-text px-0')
                    ]), className='text-center shadow-0')
                ),
                html.Button(id='button-initial-expert', className='card-button btn btn-gray m-4',
                    title='Go to the main window containing all options to configure.',
                    children=[dbc.Card(dbc.CardBody([
                        html.H5("Manual mode", className='card-title'),
                        html.Img(height='80rem', src=app.get_asset_url('icon-expert.png'),
                                 className='card-text m-3'),
                        html.Br(),
                        html.P("Fully manual configuration", className='card-text px-0')
                    ]), className='text-center shadow-0')]
                )
            ])
        ]


def initial_window_pick_band():
    """Initial window with the introduction to the EVN Observation Planner and the band selection.
    """
    return [html.Div(className='row', children=[
                html.H3("Select the observing band"),
                html.P(["At which frequency or wavelength do you want to conduct your observation? "
                    "Note that, in any case, you will still be able to change your selection "
                    "afterwards in case you want to change to or compare different bands."])
            ], style={'text:align': 'justify !important'}),
            html.Br(),
            html.Div(className='justify-content-center', children=[html.Div(
                dcc.Slider(id='initial-band', min=0, max=len(fs.bands)-1,
                       value=tuple(fs.bands).index('18cm'), step=-1,
                       marks={i: fq for i,fq in enumerate(fs.bands)},
                       persistence=True,
                       updatemode='drag', included=False)), html.Br(),
                html.Div(id='initial-pickband-label', className='row justify-content-center', children=''),
                html.Br(),
                html.Div(className='row text-center',
                         children=html.Button('Continue', id='button-pickband',
                            className='btn btn-primary btn-lg'))
            ], style={'min-width': '33rem'}),
        ]


def initial_window_pick_network(app, vlbi_networks):
    """Initial window to introduce the default VLBI network(s) to be used.
    It will only allow the ones that can observe at the given wavelength.
    """
    return [html.Div(className='row', children=[
            html.H3('Select the VLBI Network(s)'),
            html.P(["Which network(s) will be used to observe your source? ", html.Br(),
                "This will select the default antennas from each network. Note that later "
                "you will be able to add or remove antennas as wished."])
            ]),
            html.Br(),
            html.Div(className='row justify-content-center', children=[
                network_card(app, a_network, vlbi_networks[a_network]['name'], "")
                for a_network in vlbi_networks if a_network != 'e-EVN'
            ]),
            html.Br(),
            html.Div(className='row justify-content-center',
                     children=html.Button('Continue', id='button-picknetwork',
                                disabled=True, className='btn btn-primary btn-lg')),
            # html.Div(hidden=True, children=[
            #     dcc.Dropdown(id='initial-array', options=[{'label': n, 'value': n} \
            #                  for n in vlbi_networks if n != 'e-EVN'], value=[],
            #                  persistence=True, multi=True)
            # ]),
            html.Div(style={'height': '20rem'})
        ]


def initial_window_pick_time():
    """Initial (second) window with the introduction to the EVN Observation Planner and
    the option to pick a specific observing time or let the tool to find them.
    """
    return [
        html.Div(className='row justify-content-center', children=[
            html.H3('Select the observing epoch and target'),
            html.P(["Introduce the epoch when the observation will be conducted and the target "
                "source to observe, if known,, or let the tool to find "
                "the most appropriate times for a given source/network."]),
            html.Br(),
            html.Div(className='row justify-content-center', children=[dbc.FormGroup([
                dbc.RadioItems(options=[{"label": "Do not specify source nor epoch",
                                         "value": SourceEpoch.UNSPECIFIED.value},
                                        {"label": "Find best epoch for the given source",
                                         "value": SourceEpoch.ONLY_SOURCE.value},
                                        {"label": "Define observing epoch and target source",
                                         "value": SourceEpoch.SOURCE_AND_EPOCH.value}],
                               value=SourceEpoch.SOURCE_AND_EPOCH.value, id="initial-timeselection",
                               inline=True, persistence=True),
                html.Br(),
                ], className='col-6', inline=True),
                html.Div(className='col-6', children=[
                    html.Small("", id='timeselection-div-smalltext', style={'color': '#999999'})
                ])]
            ),
            html.Br(),
            html.Div(className='col-7', children=[
                html.Div(id='initial-timeselection-div-source', className='row justify-content-center',
                    hidden=True, children=[
                            html.Label(['Source (name or coordinates)',
                                *tooltip(idname='popover-target',
                                     message="Source name or coordinates. " \
                                         "You may see an error if the given name is not properly resolved. "
                                         "J2000 coordinates are assumed in both forms: 00:00:00 00:00:00 or " \
                                         "00h00m00s 00d00m00s.")
                            ]),
                            html.Div(className='form-group', children=[
                                dcc.Input(id='initial-source', value=None, type='text',
                                          className='form-control', placeholder="hh:mm:ss dd:mm:ss",
                                          persistence=True),
                                html.Small(id='initial-error_source',
                                           className='form-text text-muted'),
                            ])
                    ]),
                    html.Div(id='initial-timeselection-div-epoch', hidden=False, children=[
                        html.Label([html.Br(), 'Start of observation (UTC)']),
                        *tooltip(idname='popover-startime', message="Select the date and "
                                    "time of the start of the observation (Universal, UTC, "
                                    "time). You will also see the day of the year (DOY) in "
                                    "brackets once the date is selected."),
                        html.Br(),
                        dcc.DatePickerSingle(id='initial-starttime', date=None, min_date_allowed=dt(1900, 1, 1),
                                             max_date_allowed=dt(2100, 1, 1),
                                             display_format='DD-MM-YYYY (DDD)',
                                             placeholder='Start date',
                                             first_day_of_week=1,
                                             initial_visible_month=dt.today(),
                                             persistence=True,
                                             className='form-picker'),
                        dcc.Dropdown(id='initial-starthour', placeholder="Start time (UTC)", value=None,
                                     options=[{'label': f"{hm//60:02n}:{hm % 60:02n}", \
                                               'value': f"{hm//60:02n}:{hm % 60:02n}"} \
                                              for hm in range(0, 24*60, 15)],
                                     persistence=True, className='form-hour'),
                        html.Small(id='initial-error_starttime', style={'color': 'red'},
                                   className='form-text text-muted')
                    ]),
                    html.Div(id='initial-timeselection-div-duration', className='row justify-content-center',
                        hidden=True, children=[
                            html.Label([html.Br(), 'Duration of the observation (in hours)']),
                            html.Div(className='form-group', children=[
                                dcc.Input(id='initial-duration', value=None, type='number', className='form-control',
                                       placeholder="Duration in hours", persistence=True, inputMode='numeric'),
                                html.Small(id='initial-error_duration', className='form-text text-danger')
                            ]),
                    ]),
    #             html.Div(id='initial-timeselection-div-guess', className='row justify-content-center',
    #                 children=[
    #                     html.Small("Choose this option if you just want to find out when your source will "
    #                         "be visible. It will pick the time range when more than 3 antennas can observe "
    #                         "simultaneously."),
    #                     html.Small("Note that this option may not provide the best (expected) results in the "
    #                         "case of combining different networks that are far apart (e.g. LBA together with "
    #                         "the EVN).")
    #             ]),
    #             html.Div(id='initial-timeselection-div-epoch', children=[
    #                 html.Label('Start of observation (UTC)'),
    #                 *tooltip(idname='popover-startime', message="Select the date and "
    #                             "time of the start of the observation (Universal, UTC, "
    #                             "time). You will also see the day of the year (DOY) in "
    #                             "brackets once the date is selected."),
    #                 html.Br(),
    #                 dcc.DatePickerSingle(id='initial-starttime', date=None, min_date_allowed=dt(1900, 1, 1),
    #                                      max_date_allowed=dt(2100, 1, 1),
    #                                      display_format='DD-MM-YYYY (DDD)',
    #                                      placeholder='Start date',
    #                                      first_day_of_week=1,
    #                                      initial_visible_month=dt.today(),
    #                                      persistence=True,
    #                                      className='form-picker'),
    #                 dcc.Dropdown(id='initial-starthour', placeholder="Start time (UTC)", value=None,
    #                              options=[{'label': f"{hm//60:02n}:{hm % 60:02n}", \
    #                                        'value': f"{hm//60:02n}:{hm % 60:02n}"} \
    #                                       for hm in range(0, 24*60, 15)],
    #                              persistence=True, className='form-hour'),
    #                 html.Small(id='initial-error_starttime', style={'color': 'red'},
    #                            className='form-text text-muted'),
    #                 html.Label('Duration of the observation (in hours)'),
    #                 html.Div(className='form-group', children=[
    #                     dcc.Input(id='initial-duration', value=None, type='number', className='form-control',
    #                                placeholder="Duration in hours", persistence=True, inputMode='numeric'),
    #                     html.Small(id='initial-error_duration',
    #                                className='form-text text-danger')
    #                 ])
    #             ])
        ]),
        html.Span(style={'height': '2rem'}),
        html.Div(className='row justify-content-center',
             children=html.Button('Continue', id='button-picktimes',
                        className='btn btn-primary btn-lg')),
    ])
    ]


def initial_window_pick_mode(app):
    """Initial window to select the observing mode: spectral line or continuum observation.
    """
    return [html.Div(className='row justify-content-center', children=[
                html.H3("What type of observation do you want?"),
                html.P(["This information will only be considered to pick up the default observing "
                    "setup that is more appropriate for you (in terms of bandwidth, number of spectral "
                    "channels, etc). You will still be able to tune these parameters afterwards."])
            ], style={'text:align': 'justify !important'}),
            html.Br(),
            html.Div(className='row card-deck justify-content-center text-center', children=[
                html.Button(id='button-mode-continuum', className='card-button btn btn-gray m-4',
                    title='Helps you to configure your observation in a step-by-step process.',
                    children=dbc.Card(dbc.CardBody([
                        html.H5("Continuum mode", className='card-title'),
                        html.Img(height='80rem', src=app.get_asset_url('waves.png'),
                                 className='card-text m-3'),
                        html.Br(),
                        html.P("Provides the maximum sensitivity", className='card-text px-0'),
                        html.P("Uses the maximum bandwidth available", className='card-text px-0')
                    ]), className='text-center shadow-0')
                ),
                html.Button(id='button-mode-line', className='card-button btn btn-gray m-4',
                            title='Go to the main window containing all options to configure.',
                            children=[dbc.Card(dbc.CardBody([
                                html.H5("Spectral line mode", className='card-title'),
                                html.Img(height='80rem', src=app.get_asset_url('waves-line.png'),
                                         className='card-text m-3'),
                                html.Br(),
                                html.P("Provides higher frequency resolution", className='card-text px-0'),
                                html.P("In general uses a reduced bandwidth", className='card-text px-0')
                            ]), className='text-center shadow-0')]
                )
            ]),
            html.Div(hidden=True, children=[dbc.Checklist(id='is_line',
                                                          options=[{'label': 'line obs', 'value': False}],
                                                          value=[])])
            ]


def initial_window_final():
    """The observing mode will be either 'cont' or 'line'.
    """
    return [html.Div(children=[
        # dcc.Tabs(id='tabs', value=''),
        html.Div(className='row justify-content-center', children=[
            html.H3('You are now ready'),
            html.P(["Press compute to produce the summary for your observation. "
                    "You would then see different tabs with the information. "
                    "You will also be able to change the setup and re-compute it."]),
            html.Br(),
        ]),
        html.Div(style={'height': '4rem'}),
        html.Div(className='row justify-content-center',
                 children=html.Button('Compute', id='antenna-selection-button',
                                      className='btn btn-primary btn-lg')),
        html.Br(),
        html.Div(className='col-9 text-center justify-content-center', children=[
            dcc.Loading(id="loading", children=[html.Div(id="loading-output",
                        className='text-center justify-content-center')],
                        type="dot"),
        ]),
        html.Div(className='row text-center justify-content-center',
                 children=dcc.ConfirmDialog(id='global-error', message=''))
        ])]


