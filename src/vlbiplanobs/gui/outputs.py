# import numpy as np
# from datetime import datetime as dt
# import enum
import io
import tempfile
from pathlib import Path
from typing import Optional, Union
from datetime import datetime as dt
import numpy as np
from astropy import units as u
from astropy.time import Time
# from fpdf import FPDF
from dash import Dash, html, dcc, callback, Output, Input, State
import dash_bootstrap_components as dbc
from borb import pdf
# import plotly.express as px
from vlbiplanobs import freqsetups as fs
from vlbiplanobs import stations
from vlbiplanobs import observation
from vlbiplanobs import cli
from vlbiplanobs.gui import plots


def quantity2str(val: u.Quantity) -> str:
    return f"{val.value:.3g} {val.unit.to_string("unicode")}"


def card(children: Optional[list] = None, className: str = '', **kwargs) -> html.Div:
    return html.Div(className=' '.join(['card', className]), style={'margin': '10px'},
                    children=[html.Div(className='card-body', children=children)])


def card_result(number: str, label: str, id: str, extra_rows: Optional[list] = None,
                second_column_n: int = 0, second_column_content: Optional[list] = None):
    return card(className='bg-primary opacity-10 m-0', style={'height': '100%'}, children=[
        html.Div(className='row', children=[
            html.Div(className=f'col-{12-second_column_n} text-start text-wrap: pretty', children=[
                html.Div(className='numbers', children=[
                    html.P(className='text-white text-sm mb-0 opacity-7 font-weight-bold',
                           children=label),
                    html.H3(className='text-white font-weight-bolder mb-0 justify-content-center',
                            id=f'{id}-value', children=number)
                ])
            ]),
            html.Div(className=f'col-{second_column_n} text-end text-wrap: pretty',
                     children=second_column_content, style={'align-content': 'center'}) \
                     if second_column_n > 0 else '',
        html.Div(className='col-12', children=extra_rows) if extra_rows is not None else ''])
    ])


def download_pdf(app, pdf_filename: str) -> html.Div:
    """Buttom to download the PDF report
    """
    html.Div(html.A('Download summary as PDF', id='request-printable-version',
                    download='evn_planobs_summary.pdf', href=app.get_asset_url(pdf_filename),
                    className='btn btn-link btn-lg', style={'text-align': 'center'}),
            className='col-12 justify-content-center', style={'text-align': 'center'})


def message_card(title: str, body: str, mode: str, icon: Optional[str] = None):
    """Shows a card with a title and a body, with an icon and a mode that defines the color of the card:
    'info', 'warning', 'danger'. If 'icon' is not provided, it will take an icon by default depending on the mode
    """
    assert mode in ('danger', 'info', 'warning')
    if icon is None:
        match mode:
            case 'danger':
                icon = 'fa-triangle-exclamation'
            case 'info':
                icon = 'fa-info'
            case 'warning':
                icon = 'fa-exclamation'

    colors = {'danger': '#ef4444', 'warning': '#eab308', 'info': '#A6D4E8'}

    return html.Div(className=f'card col-12 px-4 pb-0 m-2 bg-{mode} opacity-10 border-radius-lg',
            children=html.Div(className='card-body p-2 position-relative',
                children=[html.Div(className='row align-items-center', children=[
                          html.Div(className='col-2 pl-3 text-center',
                                   children=html.Div(html.I(className=f'fa-solid {icon} '
                                                                     'text-light '
                                                                     'text-xl opacity-40',
                                                            style={'font-size': '3rem',
                                                                   'color': colors[mode]}),
                                                     )),
                          html.Div(className='col-10 text-start align-items-center',
                                   children=[html.H5(className='text-white font-weight-bolder mb-0 mt-0',
                                                     children=title),
                                             html.Span(className='text-white text-sm',
                                                       children=body)])
                   ])
                ]))


def warning_card():
    # TODO: reduce the code un the sun, phase referencing cards to this one.
    pass


def error_card(title: str, body: str) -> html.Div:
    return message_card(title, body, 'danger')


def info_card(title: str, body: str) -> html.Div:
    return message_card(title, body, 'info')


def warning_phase_referencing_high_freq(o: Optional[cli.VLBIObs] = None) -> html.Div:
    """For observations at high frequencies, it shows a warning mentioning that
    phase referencing is not possible.
    """
    if o is None:
        return html.Div()

    if o.frequency < 70*u.GHz:
        return html.Div()
    else:
        return html.Div(className='card pb-0 m-2 bg-warning opacity-10 border-radius-lg',
            children=html.Div(className='card-body p-2 position-relative',
                children=[html.Div(className='row align-items-center', children=[
                          html.Div(className='col-2 pl-3 text-center',
                                   children=html.Div(html.I(className='fa-solid fa-repeat '
                                                                     'text-light '
                                                                     'text-xl opacity-40',
                                                            style={'font-size': '3rem',
                                                                   'color': '#eab308'}),
                                                     )),
                                   # children=html.Div(className='icon icon-shape bg-white shadow '
                                   #                   'text-center border-radius-2xl',
                                   #                   children=html.I(className='fa-solid fa-repeat '
                                   #                                   'text-dark '
                                   #                                   'text-gradient text-xl opacity-10'))),
                          html.Div(className='col-10 text-start align-items-center',
                                   children=[html.H5(className='text-white font-weight-bolder mb-0 mt-0',
                                                     children="Note that phase referencing is not "
                                                              "possible"),
                                             html.Span(className='text-white text-sm',
                                                       children="At these high frequencies, your target "
                                                                "must be bright enough to directly fringe "
                                                                "on it.")])
                   ])
                ]))


def summary_freq_res(o: Optional[cli.VLBIObs] = None) -> html.Div:
    if o is None:
        return html.Div()

    try:
        chan_f = cli.optimal_units(o.bandwidth/(o.subbands*o.channels),
                          [u.GHz, u.MHz, u.kHz, u.Hz])
        vel = cli.optimal_units((2.9979e5*u.km/u.s)*(o.bandwidth/(2*o.frequency)).decompose(),
                            [u.km/u.s, u.m/u.s])
        chan_v = cli.optimal_units(2*vel/(o.subbands*o.channels),
                          [u.km/u.s, u.m/u.s])
    except:
        return html.Div()

    return card_result(f"{chan_f.value:.3g} {chan_f.unit.to_string('unicode')}",
                       f"Frequency resolution", id='freq-res-frequency',
                       second_column_n=6, second_column_content=[
                            html.Div(className='numbers', children=[
                                html.P(className='text-white text-sm mb-0 opacity-7 font-weight-bold',
                                       children="Velocity resolution"),
                                html.H3(className='text-white font-weight-bolder mb-0 justify-content-center',
                                        id='freq-res-velocity',
                                        children=f"{chan_v.value:.3g} "
                                                 f"{chan_v.unit.to_string('unicode')}")
                            ])
                       ], extra_rows=[
                            html.Br(),
                            html.Div(className='col-12', style={'color': 'var(--bs-gray-100)'},
                                     children=html.P(f"The total bandwidth of {o.bandwidth} is divided in "
                                                     f"{o.subbands} x {o.bandwidth/o.subbands} subbands, "
                                                     f"with {o.channels} spectral channels each.") \
                                              if o.subbands > 1 else \
                                              html.P(f"The total bandwidth of {o.bandwidth} is recorded in "
                                                     "a single subband, with {o.channels} spectral channels."))
                       ])


def field_of_view(o: Optional[cli.VLBIObs] = None) -> html.Div:
    if o is None:
        return html.Div()

    if o.bandwidth_smearing() is None or o.time_smearing() is None:
        return html.Div()

    bw_smearing = cli.optimal_units(o.bandwidth_smearing(), [u.deg, u.arcmin, u.arcsec])
    tm_smearing = cli.optimal_units(o.time_smearing(), [u.deg, u.arcmin, u.arcsec])
    return card_result(f"{bw_smearing.value:.2n} {bw_smearing.unit.to_string('unicode')}",
                       f"Bandwidth smearing", id='fov',
                       second_column_n=6, second_column_content=[
                            html.Div(className='numbers', children=[
                                html.P(className='text-white text-sm mb-0 opacity-7 font-weight-bold',
                                       children="Time smearing"),
                                html.H3(className='text-white font-weight-bolder mb-0 justify-content-center',
                                        id=f'{id}-value',
                                        children=f"{tm_smearing.value:.2n} "
                                                 f"{tm_smearing.unit.to_string('unicode')}")
                            ])
                       ], extra_rows=[
                            html.Br(),
                            html.Div(className='col-12', style={'color': 'var(--bs-gray-100)'},
                                     children=html.P("The field of view (FoV) will be limited to this radius"
                                                     " after correlation due to time and frequency smearing, "
                                                     "considering a 10% loss."))
                       ])

def rms(o: Optional[cli.VLBIObs] = None) -> html.Div:
    if o is None:
        return html.Div(dbc.Modal(id="sensitivity-baseline-modal", is_open=False,
                                  children=[html.Div(html.Button('View sensitivity per baseline',
                                                         id='button-sensitivity-baseline'))
                                  ]))

    try:
        if isinstance(thermal_noise := o.thermal_noise(), dict):
            rms = list(thermal_noise.values())[0]
        else:
            rms = thermal_noise

        out_rms = [quantity2str(cli.optimal_units(rms,
                                                  [u.MJy/u.beam, u.kJy/u.beam, u.Jy/u.beam,
                                                   u.mJy/u.beam, u.uJy/u.beam])),
                   quantity2str(cli.optimal_units(rms/np.sqrt(1*u.min/ \
                                                  (o.duration \
                                                  if o.duration is not None else 24*u.h)),
                                                  [u.MJy/u.beam, u.kJy/u.beam, u.Jy/u.beam,
                                                   u.mJy/u.beam, u.uJy/u.beam])),
                   quantity2str(cli.optimal_units(rms*np.sqrt(o.subbands* \
                                                  o.channels),
                                                  [u.MJy/u.beam, u.kJy/u.beam, u.Jy/u.beam,
                                                   u.mJy/u.beam, u.uJy/u.beam])),
                   show_baseline_sensitivities(o)]
    except:
        return html.Div(dbc.Modal(id="sensitivity-baseline-modal", is_open=False,
                                  children=[html.Div(html.Button('View sensitivity per baseline',
                                                         id='button-sensitivity-baseline'))
                                  ]))

    src = list(o.ontarget_time.keys())[0]
    return card_result(out_rms[0], f'rms thermal noise (for {o.ontarget_time[src]:.2g} on-target)', id='rms',
                       extra_rows=[html.Br(), html.Div(className='row', children=[
                            html.Div(className='col-6 text-start px-0 pb-2', children=[
                                html.H5(className='mb-0 font-weight-bolder text-light',
                                       id='rms-per-channel-value',
                                       children=out_rms[2]),
                                html.Label(className='text-sm mb-0 text-light',
                                       children="per spectral channel", id='tooltip-rms-channel'),
                                ], style={'text-wrap': 'pretty'}),
                                dbc.Tooltip("This is the theoretical rms noise to obtain under ideal conditions "
                                            "in a single spectral channel when integrating over the whole "
                                            " observation.", target='tooltip-rms-channel'),
                            html.Div(className='col-6 text-end px-0', children=[
                                html.H5(className='mb-0 font-weight-bolder text-light',
                                       id='rms-per-time-value',
                                       children=out_rms[1]),
                                html.Label(className='text-sm mb-0 text-light',
                                       children="on 1-min integration", id='tooltip-rms-time'),
                                dbc.Tooltip("This is the theoretical rms noise to obtain under ideal conditions "
                                            "when integrating the data from the full bandwidth with one minute "
                                            "time integration.", target='tooltip-rms-time'),
                            ], style={'text-wrap': 'pretty'}),
                                ]),
                                html.Div(className='row', children=[
                                    html.Div(html.Button('View sensitivity per baseline',
                                                         id='button-sensitivity-baseline',
                                                         # color='light',
                                                         # className='btn bg-light btn-lg ' \
                                                         #           'mx-auto w-75 m-4 p-2 mb-0 active',
                                                         className='btn btn-white btn-sm w-100 '
                                                         'mb-0 active',
                                                         style={'position': 'sticky', 'top': '20px',
                                                                'color': '#9DB7C4', 'box-shadow': 'none'}),
                                        )
                                   ]),
                                dbc.Modal(id='sensitivity-baseline-modal',
                                          size='xl', is_open=False,
                                          children=show_baseline_sensitivities(o)),
                                # html.Div(className='col-12', id='table-sensitivities',
                                #          children=None)
                        ])


def show_baseline_sensitivities(o: Optional[cli.VLBIObs] = None) -> html.Div:
    if o is None:
        return html.Div([dbc.ModalHeader(dbc.ModalTitle("Sensitivity per baseline")),
                          dbc.ModalBody("")], id='sens-baseline-style', style={'display': 'none'})

    table_header = [html.Thead(html.Tr([html.Th("")] + [html.Th(s.codename,
                                                        className='text-left text-sm p-0 pb-2 ps-3') \
                                                        for s in o.stations]))]
    all_baselines = np.unique(np.array([v.to(u.mJy/u.beam).value \
                    for v in o.baseline_sensitivity().values()]))
    lower_div, higher_div = np.quantile(all_baselines, 0.33), np.quantile(all_baselines, 0.66)
    rows = []
    for i, s in enumerate(o.stations):
        # a_row = [html.Td(" ")]*i
        a_row = [html.Td(" ")]*i
        for j in range(i, len(o.stations)):
            val = o.baseline_sensitivity(s.codename,
                                         o.stations[j].codename).to(u.mJy/u.beam).value
            val_badge = 'bg-success' if val < lower_div else 'bg-danger' \
                        if val > higher_div else 'bg-warning'
            a_row.append(html.Td(html.Span(f"{val:4.2f}", className='badge ' + val_badge)))

        rows.append(html.Tr([html.Th(s.codename, className='text-left text-sm p-0 pb-2 ps-3')] + \
            a_row))

    table = dbc.Table(table_header + [html.Tbody(rows)], bordered=True,
                      className='table align-items-center mb-0')
    return html.Div([dbc.ModalHeader(dbc.ModalTitle("Sensitivity per baseline "
                                           f"({(u.mJy/u.beam).to_string("unicode")}) for one-minute "
                                           "time integration")),
                     html.P("The following table shows the sensitivity for each baseline (or auto-corrleation) "
                            "for one-minute time integration (considering the full bandwidth available for "
                            "each baseline. The colors highlight the most sensitive baselines (green color) "
                            "to the less sensitive baselines (red colors).", className='text-dark p-2'),
                     dbc.ModalBody(table, className='table-responsive')], id='sens-baseline-style',
                    style={'display': 'block'})


def resolution(o: Optional[cli.VLBIObs] = None) -> html.Div:
    if o is None:
        return html.Div()

    if not o.sourcenames:
        synth_beam = o.synthesized_beam()[list(o.synthesized_beam().keys())[0]]
    else:
        synth_beam = o.synthesized_beam()[o.sourcenames[0]]

    bmaj = cli.optimal_units(synth_beam['bmaj'], [u.deg, u.arcmin, u.arcsec, u.mas, u.uas])
    bmin = synth_beam['bmin'].to(bmaj.unit)
    bmin_elip = max(int(bmin.value*80/bmaj.value), 2)
    return card_result([f"{bmaj.value:2.1f} x {quantity2str(bmin)}", html.Sup("2"),
                        f", {synth_beam['pa'].value:2.0f}º"], 'Angular Resolution', id='res',
                       # second_column_n=3, second_column_content=[
                       # ellipse(bmaj="32px", bmin="10px", pa="10deg")
                       # TODO: double check angle. Is it x-inverted??????
                       extra_rows=[html.Br(), html.Div(ellipse(bmaj="80px",
                                                               bmin=f"{bmin_elip}px",
                                                               pa=f"{synth_beam['pa'].value}deg"),
                                                       className='row justift-content-center py-4',
                                                       style={'display': 'flex',
                                                              'justify-content': 'center'})])

def spectral_resolution(o: Optional[cli.VLBIObs] = None) -> html.Div:
    if o is None:
        return html.Div()

    vel = cli.optimal_units((2.9979e5*u.km/u.s)*(o.bandwidth/ \
                            (o.frequency*o.subbands*o.channels)).decompose(),
                            [u.km/u.s, u.m/u.s])
    chan_width = cli.optimal_units(o.bandwidth/(o.subbands*o.channels).decompose(),
                                   [u.GHz, u.MHz, u.kHz, u.Hz])

    return card_result(quantity2str(vel), 'Spectral Resolution', id='res',
                       extra_rows=[html.Div(className='row', children=[
                        html.Div(className='col-12', children=[
                            html.Label(f"Per spectral channel ({quantity2str(chan_width)} wide)"),
                            html.Div(className='numbers', children=[
                                html.P(className='text-white text-sm mb-0 opacity-7 font-weight-bold',
                                       children='Velocity Resolution'),
                                html.H3(className='text-white font-weight-bolder mb-0 ' \
                                                  'justify-content-center',
                                        id='vel-res-value', children=quantity2str(vel))
                            ])
                        ])
                       ])])


def ellipse(bmaj, bmin, pa, color='white', z_index=1, position='relative', margin_top='', className=''):
    """Returns a html.Div element that draws an ellipse with a semimajor axis bmaj,
    semiminor axis bmin, and position angle (as defined in radio astronomy) pa.
    bmaj,bmin, pa must be strings recognized by HTML/CSS.
    """
    return html.Div(children=[], style={'width': bmaj, 'height': bmin,
                    'border-radius': '50%', 'background': color, 'position': position,
                    'transform': f"rotate({pa})", 'z-index': z_index, 'justify-content': 'center',
                    'vertical-align': 'middle', 'margin-top': margin_top}, className=className)


def card_sun_warning(head_line: str, body_line: str, hidden: bool = False) -> html.Div:
    return html.Div(className='card pb-0 m-2 bg-warning opacity-10 border-radius-lg', hidden=hidden,
        children=html.Div(className='card-body p-2 position-relative',
            children=[html.Div(className='row align-items-center', children=[
                      html.Div(className='col-2 pl-3 text-center',
                               children=html.Div(html.I(className='fa-solid fa-sun '
                                                                  'text-light '
                                                                  'text-xl opacity-40',
                                                        style={'font-size': '3em', 'color': '#eab308'}))),
                      html.Div(className='col-10 text-start align-items-center',
                               children=[html.H5(className='text-white font-weight-bolder mb-0 mt-0',
                                                 children=head_line),
                                         html.Span(className='text-white text-sm',
                                                   children=body_line)])
               ])
            ]))


def card_ant_warning(head_line: str, body_line: str, hidden: bool = False) -> html.Div:
    return html.Div(className='card pb-0 m-2 bg-warning opacity-10 border-radius-lg', hidden=hidden,
        children=html.Div(className='card-body p-2 position-relative',
            children=[html.Div(className='row align-items-center', children=[
                      html.Div(className='col-2 pl-3 text-center',
                               children=html.Div(className='icon icon-shape bg-white shadow '
                                                 'text-center border-radius-2xl',
                                                 children=html.I(className='fa-solid fa-satellite-dish '
                                                                 'text-dark '
                                                                 'text-xl opacity-10'))),
                      html.Div(className='col-10 text-start align-items-center',
                               children=[html.H5(className='text-white font-weight-bolder mb-0 mt-0',
                                                 children=head_line),
                                         html.Span(className='text-white text-sm',
                                                   children=body_line)])
               ])
            ]))


def ant_warning(o: Optional[cli.VLBIObs]) -> html.Div:
    """Returns a warning card if some antennas selected in the observation cannot observe.
    """
    if o is None:
        return html.Div()

    dropped_stations = [s for s in o.stations if o.band not in s.bands]
    if len(dropped_stations) > 0:
        ("[yellow]The following antennas were ignored"
               f"because they cannot observe at {band}: {', '.join(dropped_stations)}[/yellow]")

        return card_ant_warning("Some antennas cannot observe at this band",
                                f"The antenna{'s' if len(o.antennas) > 1 else ''} "
                                f"{', '.join([a.name for a in o.antennas])} cannot observe at {band}.")

    return html.Div()


def sun_warning(o: Optional[cli.VLBIObs] = None) -> html.Div:
    """Returns a warning card if the Sun is too close to the target source.
    Otherwise it returns an empty Div.
    """
    if o is None:
        return html.Div()

    sun_const = o.sun_constraint()
    sun_limit = o.sun_limiting_epochs()
    if o.times is None:
        if len(sun_limits := list(sun_limit.values())[0]) > 0:
            sun_const = sun_const[list(sun_const.keys())[0]]
            t0, t1 = sun_limits[0].datetime, sun_limits[-1].datetime
            if t0 == t1:
                return card_sun_warning("The Sun gets too close to the source",
                                        f"On {t0.strftime('%d %b')} "
                                        f"(with a minimum separation of {sun_const.value:.0f}"
                                        f"{sun_const.unit.to_string('unicode')}).")
            elif t0.month == t1.month:
                return card_sun_warning("The Sun gets too close to the source",
                                        f"During {t0.day}-{t1.day} {t0.strftime('%B')} "
                                        f"(with a minimum separation of {sun_const.value:.0f}"
                                        f"{sun_const.unit.to_string('unicode')}).")
            elif t0.year == t1.year:
                return card_sun_warning("The Sun gets too close to the source",
                                        f"From {t0.strftime('%d %B')} to {t1.strftime('%d %B')} "
                                        f"(with a minimum separation of {sun_const.value:.0f}"
                                        f"{sun_const.unit.to_string('unicode')}).")
            else:
                return card_sun_warning("The Sun gets too close to the source",
                                        f"From {t0.strftime('%d %B %Y')} to {t1.strftime('%d %B %Y')} "
                                        f"(with a minimum separation of {sun_const.value:.0f}"
                                        f"{sun_const.unit.to_string('unicode')}).")
    else:
        if len(sun_separation := list(sun_const.values())) > 0 and sun_separation[0] is not None:
            sun_const = sun_const[list(sun_const.keys())[0]]
            return card_sun_warning("The Sun is too close to the source!",
                                    f"With a minimum separation of {sun_separation[0].value:.0f}"
                                    f"{sun_const.unit.to_string('unicode')} during the observation.")

    return html.Div()


def plot_elevations(o: Optional[cli.VLBIObs] = None) -> html.Div:
    if o is None:
        return html.Div()

    return card([html.Div(className='card-header pb-0', children=html.H5('Source Elevation')),
                 dcc.Graph(id='fig-elevations', figure=plots.elevation_plot(o)),
                 html.Div(print_observability_ranges(o))],
                className='')


def print_observability_ranges(o: cli.VLBIObs) -> html.Div:
    if o is None:
        return html.Div()

    text = []
    if o.times is None:
        srcup = o.is_observable_at(Time('2025-09-21', scale='utc') + np.arange(0.0, 1.005, 0.01)*u.day)
    else:
        srcup = o.is_observable()

    ablockname = list(srcup.keys())[0]
    text += [html.Br()]
    if o.times is None:
        when_everyone = o.when_is_observable(mandatory_stations='all',
                                             return_gst=True)[ablockname]
        if len(when_everyone) > 0:
            text += ["Everyone can observe the source at: " + \
                    ', '.join([t1.to_string(sep=':', fields=2, pad=True) + '-' +
                    t2.to_string(sep=':', fields=2, pad=True) +
                    ' GST.' for t1, t2 in when_everyone])]
        else:
            text += ["The source cannot be observed by all stations at the same time."]
    else:
        srcupalways = o.is_always_observable()
        when_everyone = o.when_is_observable(mandatory_stations='all')[ablockname]
        if len(when_everyone) > 0:
            text += ["Everyone can observe the source on " + \
                    ', '.join([t1.strftime('%d %b %Y %H:%M')+'-'+t2.strftime('%H:%M') +
                    ' UTC.' for t1, t2 in when_everyone])]
        else:
            text += ["The source cannot be observed by all stations at the same time."]

        text += [html.Br()]
        if any(srcupalways[ablockname].values()):
            ant_can = [ant for ant, b in srcupalways[ablockname].items() if b]
            if len(ant_can) <= len(o.stations)/2:
                text += [f"Only {', '.join(ant_can)} can observe the source at all times."]
            else:
                text += ["All antennas can observe the source at all times except " + \
                        ', '.join([ant for ant in o.stations.station_codenames if ant not in ant_can]) + '.']
        else:
            text += ["But there are no antennas that can observe the source at all time."]

    text += [html.Br()]
    min_stat = 3 if len(o.stations) > 3 else min(2, len(o.stations))
    if len(o.stations) > 2:
        text += [f"The optimal visibility range (> {min_stat} antennas) occurs on "]
        if o.times is None:
            text += [', '.join([t1.to_string(sep=':', fields=2, pad=True) + '-' +
                               (t2 + (24*u.hourangle if np.abs(t1 - t2) < 0.1*u.hourangle
                                      else 0.0*u.hourangle)).to_string(sep=':', fields=2, pad=True) +
                               ' GST.' for t1, t2
                               in o.when_is_observable(min_stations=min_stat,
                                                       return_gst=True)[ablockname]])]
        else:
            text += [', '.join([t1.strftime('%d %b %Y %H:%M')+'-'+t2.strftime('%H:%M') +
                               ' UTC.' for t1, t2
                               in o.when_is_observable(min_stations=min_stat)[ablockname]])]

    return text


def plot_uv_coverage(o: Optional[cli.VLBIObs] = None) -> html.Div:
    if o is None:
        return html.Div()

    longest_bl = o.longest_baseline()[list(o.longest_baseline())[0]]
    ant_l1, ant_l2 = longest_bl[0].split('-')
    # Using dummy units to allow the conversion
    longest_bl_lambda = longest_bl[1]/o.wavelength
    longest_bl_lambda = cli.optimal_units(longest_bl_lambda*u.m, [u.Gm, u.Mm, u.km])
    shortest_bl = o.shortest_baseline()[list(o.longest_baseline())[0]]
    ant_s1, ant_s2 = shortest_bl[0].split('-')
    # Using dummy units to allow the conversion
    shortest_bl_lambda = shortest_bl[1]/o.wavelength
    shortest_bl_lambda = cli.optimal_units(shortest_bl_lambda*u.m, [u.Gm, u.Mm, u.km])

    return card([html.Div(className='card-header pb-0', children=html.H5('(u, v) Coverage')),
                 dcc.Graph(id='fig-coverage', figure=plots.uvplot(o)),
                 html.P(""),
                 html.Label("Highlight antennas:", id='select-ant-uv-label',
                            htmlFor='select-ant-uv-plot'),
                 dcc.Dropdown(multi=True, id="select-antenna-uv-plot",
                              options=[{'label': f"{ant.name} ({ant.codename})", 'value': ant.codename} \
                                       for ant in o.stations]),
                 html.Br(), html.Br(),
                 html.Div(className='row', children=[
                     html.Div(className='col-6 text-start text-wrap: pretty', children=[
                        html.Label(className='text-sm mb-0',
                               children="Longest baseline"),
                        html.H5(className='mb-0 font-weight-bolder',
                               id='baseline-long-value',
                               children=f"{cli.optimal_units(longest_bl[1], [u.km, u.m]):.5n} "
                                        f"({longest_bl_lambda.value:.3n} "
                                        f"{longest_bl_lambda.unit.name[0]}λ)"),
                        html.Label(className='text-sm mb-0',
                               children=f"{o.stations[ant_l1].name} - "
                                        f"{o.stations[ant_l2].name} ({ant_l1}-{ant_l2})"),
                     ]),
                     html.Div(className='col-6 text-end text-wrap: pretty', children=[
                        html.Label(className='text-sm mb-0',
                               children="Shortest baseline"),
                        html.H5(className='mb-0 font-weight-bolder',
                               id='baseline-short-value',
                               children=f"{cli.optimal_units(shortest_bl[1], [u.km, u.m]):.5n} "
                                        f"({shortest_bl_lambda.value:.3n} "
                                        f"{shortest_bl_lambda.unit.name[0]}λ)"),
                        html.Label(className='text-sm mb-0',
                               children=f"{o.stations[ant_s1].name} - "
                                        f"{o.stations[ant_s2].name} ({ant_s1}-{ant_s2})"),
                     ]),

                 ])
                 ], className='')


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

def baselines(app, o: cli.VLBIObs) -> html.Div:
    """Returns the card showing the longest and shortest baseline in the array
    """
    return card([html.Div(className='card-header pb-0', children=html.H5('Baseline')),
                 html.Div()],
                className='')




# THIS IS A GOOD REFERENCE FOR THE DOCUMENTATION LINK
# <div class="card card-background shadow-none card-background-mask-secondary" id="sidenavCard">
#         <div class="full-background" style="background-image: url('../assets/img/curved-images/white-curved.jpg')"></div>
#         <div class="card-body text-start p-3 w-100">
#           <div class="icon icon-shape icon-sm bg-white shadow text-center mb-3 d-flex align-items-center justify-content-center border-radius-md">
#             <i class="ni ni-diamond text-dark text-gradient text-lg top-0" aria-hidden="true" id="sidenavCardIcon"></i>
#           </div>
#           <div class="docs-info">
#             <h6 class="text-white up mb-0">Need help?</h6>
#             <p class="text-xs font-weight-bold">Please check our docs</p>
#             <a href="https://www.creative-tim.com/learning-lab/bootstrap/license/soft-ui-dashboard" target="_blank" class="btn btn-white btn-sm w-100 mb-0">Documentation</a>


def worldmap_plot(o: Optional[cli.VLBIObs] = None) -> html.Div:
    if o is None:
        return html.Div()

    data = {"lat": [], "lon": [], "color": [], "symbol": [], "mode": "markers",
            "name": [], "text": [], "hovertemplate": [], "observes": []}
    try:
        ant_observes = o.can_be_observed()[list(o.can_be_observed().keys())[0]]
    except (ValueError, IndexError):
        ant_observes = {ant.codename: True for ant in o.stations}

    for ant in o.stations:
        data["lat"].append(ant.location.lat.value)
        data["lon"].append(ant.location.lon.value)
        data["name"].append(ant.name)
        data["observes"].append(ant_observes[ant.codename])
        data["text"].append(f"{ant.name}<br>({ant.country})<br> {ant.diameter}")
        data["hovertemplate"].append(f"{ant.name}<br>({ant.country})<br> {ant.diameter}<extra></extra>")

    return card(html.Div(className='row justify-content-center',
                    children=[dcc.Graph(figure=plots.plot_worldmap_stations(data),
        responsive=True, className='fig-on-card',
                  config={'frameMargins': 0, 'showLink': False, 'displaylogo': False})]),
                className='col-12 mx-0 px-0')
                  # config={'frameMargins': -10, 'showLink': False, 'displaylogo': False}), html.Br()]))

def button_summary(o: cli.VLBIObs) -> html.Div:
    return dbc.Button(className='btn-lg mx-auto w-75 m-4 p-2 btn btn-outline-secondary',
                      color='secondary', outline=True,
                                id='button-download', style={'position': 'sticky', 'top': '20px'},
                                children="DOWNLOAD SUMMARY AS PDF", download='planobs_summary.pdf')


def summary_pdf(o: cli.VLBIObs):
    """Creates a PDF file with the summary of the observation and includes the elevation plot figure.
    """
    buffer = io.BytesIO()
    doc: pdf.Document = pdf.Document()
    page = pdf.Page()
    doc.add_page(page)
    layout: pdf.PageLayout = pdf.SingleColumnLayout(page)
    layout.add(pdf.Paragraph("EVN Observation Planner - Summary Report", font_size=20, font='Helvetica-bold',
                         horizontal_alignment=pdf.Alignment.CENTERED))
    text = f"Observation to be conducted at {o.band.replace('cm', ' cm')}"
    if o.times is not None:
        if o.times[0].datetime.date() == o.times[-1].datetime.date():
            layout.add(pdf.Paragraph(f"{text} from {o.times[0].strftime('%d %b %Y %H:%M')}–"
                                     f"{o.times[-1].strftime('%H:%M')} UTC."))
        elif (o.times[-1] - o.times[0]) < 24*u.h:
            layout.add(pdf.Paragraph(f"{text} from {o.times[0].strftime('%d %b %Y %H:%M')}–"
                                     f"{o.times[-1].strftime('%H:%M')} (+1d) UTC."))
        else:
            layout.add(pdf.Paragraph(f"{text} from {o.times[0].strftime('%d %b %Y %H:%M')} to "
                                     f"{o.times[-1].strftime('%d %b %H:%M')} UTC."))
    else:
        min_stat = 3 if len(o.stations) > 3 else min(2, len(o.stations))
        if len(o.stations) > 2:
            if o.times is not None:
                srcup = o.is_observable()
            else:
                srcup = o.is_observable_at(Time('2025-09-21', scale='utc') + np.arange(0.0, 1.005, 0.01)*u.day)

            for ablockname, antbool in srcup.items():
                gst_range = (', '.join([t1.to_string(sep=':', fields=2, pad=True) + '--' +
                                  (t2 + (24*u.hourangle if np.abs(t1 - t2) < 0.1*u.hourangle
                                         else 0.0*u.hourangle)).to_string(sep=':', fields=2, pad=True) +
                                  ' GST.' for t1, t2
                                  in o.when_is_observable(min_stations=min_stat,
                                                          return_gst=True)[ablockname]]))
                layout.add(pdf.Paragraph(f"{text}. Optimal visibility range (> {min_stat} antennas) at "
                                         f"{gst_range}{' for '+ablockname if len(srcup) > 1 else ''}"))

    sun_const = o.sun_constraint()
    sun_limit = o.sun_limiting_epochs()
    for ablockname in o.sun_limiting_epochs():
        if o.times is None:
            if len(sun_limit[ablockname]) > 0:
                text = "Note the the Sun is too close to this source"
                t0, t1 = sun_limit[ablockname][0].datetime, sun_limit[ablockname][-1].datetime
                if t0 == t1:
                    text += f" on {t0.strftime('%d %b')}"
                elif t0.month == t1.month:
                    text += f" on {t0.day}-{t1.day} {t0.strftime('%b')}"
                elif t0.year == t1.year:
                    text += f" from {t0.strftime('%d %b')} to {t1.strftime('%d %b')}"
                else:
                    text += f" from {t0.strftime('%d %b %Y')} to {t1.strftime('%d %b %Y')}"
                text += f" (minimum separation of {sun_const[ablockname]:.02f})"
                text += f" for {ablockname}." if len(o.sun_limiting_epochs()) > 1 else "."
                layout.add(pdf.Paragraph(text, font_color=pdf.HexColor("#FF0000")))
        else:
            if sun_const[ablockname] is not None:
                text = "Note the the Sun is too close to the source"
                text += f" {ablockname}" if len(o.sun_limiting_epochs()) > 1 else " "
                text += f"({sun_const[ablockname]:.02f} away)."
                layout.add(pdf.Paragraph(text, font_color=pdf.HexColor("#FF0000")))

    if o.duration is not None:
        layout.add(pdf.Paragraph("With a total duration of "
                                 f"{cli.optimal_units(o.duration, [u.h, u.min, u.s]):.01f} "
                                 f"({cli.optimal_units(o.ontarget_time[list(o.ontarget_time.keys())[0]],
                                                  [u.h, u.min, u.s]):.01f} on target). "
                                 f"Total output FITS file size: {o.datasize():.02f}."))

    layout.add(pdf.Paragraph(f"Participating stations: {', '.join(o.stations.station_codenames)}."))
    if len(o.scans) > 0:
        for ablock in o.scans.values():
            layout.add(pdf.Paragraph(f"Target source: {'\n'.join([s.name + ' (' + s.coord.to_string('hmsdms') \
                                                              + ').' for s in ablock.sources()])}"))
    else:
        layout.add(pdf.Paragraph("No sources defined."))

    if None not in (o.datarate, o.bandwidth, o.subbands):
        val = cli.optimal_units(o.datarate, [u.Gbit/u.s, u.Mbit/u.s])
        layout.add(pdf.Paragraph(f"\nData rate of {val:.0f}, "
               f"producing a total bandwidth of {cli.optimal_units(o.bandwidth, [u.MHz, u.GHz])}, "
               f" divided in {o.subbands} x {int(o.bandwidth.value/o.subbands)}-"
               f"{o.bandwidth.unit} subbands, with {o.channels} channels each, "
               f"{o.polarizations} polarization, and {o.inttime:.01f} integration time."))
    else:
        layout.add(pdf.Paragraph("No setup (data rate, bandwidth, number of subbands) specified."))

    for ablock in o.scans.values():
        for src in o.sources():
            if len(o.scans) > 1:
                layout.add(pdf.Paragraph(f"For the source {src.name}", font='Helvetica-bold'))

            rms = cli.optimal_units(o.thermal_noise()[src.name], [u.Jy/u.beam, u.mJy/u.beam, u.uJy/u.beam])
            rms_chan = cli.optimal_units(rms/np.sqrt(1*u.min/ \
                                         (o.duration if o.duration is not None else 24*u.h)),
                                         [u.MJy/u.beam, u.kJy/u.beam, u.Jy/u.beam,
                                          u.mJy/u.beam, u.uJy/u.beam])
            rms_min = cli.optimal_units(rms*np.sqrt(o.subbands*o.channels),
                                        [u.MJy/u.beam, u.kJy/u.beam, u.Jy/u.beam,
                                         u.mJy/u.beam, u.uJy/u.beam])
            layout.add(pdf.Paragraph(f"Thermal rms noise: "
                                     f"{rms:.3g} ({rms_chan:.3g} per spectral "
                                     f"channel and {rms_min:.3g} per one-minute time integration."))
            synth_beam = o.synthesized_beam()[src.name]
            bmaj = cli.optimal_units(synth_beam['bmaj'], [u.deg, u.arcmin, u.arcsec, u.mas, u.uas])
            bmin = synth_beam['bmin'].to(bmaj.unit)
            bmin_elip = max(int(bmin.value*80/bmaj.value), 2)
            layout.add(pdf.Paragraph(f"Synthesized beam{' for ' + src.name if len(o.scans.values()) \
                                                        > 1 else ''}: "
                                     f"{bmaj.value:2.1f} x {bmin:2.1f}"
                                     f", {synth_beam['pa'].value:2.0f}º."))

    bw_smearing = cli.optimal_units(o.bandwidth_smearing(), [u.deg, u.arcmin, u.arcsec])
    tm_smearing = cli.optimal_units(o.time_smearing(), [u.deg, u.arcmin, u.arcsec])
    layout.add(pdf.Paragraph(f"Field of view limited to {bw_smearing:.2g} (from frequency smearing) "
                             f"and {tm_smearing:.2g} (from time smearing), considering 10% loss."))

    if len(o.scans.keys()) > 0:
        fig = plots.elevation_plot(o, show_colorbar=True)
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tempfig:
            fig.write_image(tempfig.name)
            figpath = Path(tempfig.name)

        layout.add(pdf.Image(figpath, width=414, height=265, horizontal_alignment=pdf.Alignment.CENTERED))

    pdf.PDF.dumps(buffer, doc)
    buffer.seek(0)
    return buffer

            # if self.times is not None or self.duration is not None:
            #     rprint("\n[bold green]Expected outcome[/bold green]:")
            #     rprint("[dim](for a +/- 45° elevation source)[/dim]")


    def get_fig_dirty_map(self):
        raise NotImplementedError
        # Right now I only leave the natural weighting map (the uniform does not always correspond to the true one)
        dirty_map_nat, laxis = self.get_dirtymap(pixsize=1024, robust='natural', oversampling=4)
        # TODO: uncomment these two lines and move them outside observation. Flexibility
        # fig1 = px.imshow(img=dirty_map_nat, x=laxis, y=laxis[::-1], labels={'x': 'RA (mas)', 'y': 'Dec (mas)'}, \
        #         aspect='equal')
        # fig = make_subplots(rows=1, cols=1, subplot_titles=('Natural weighting',), shared_xaxes=True, shared_yaxes=True)
        fig.add_trace(fig1.data[0], row=1, col=1)
        mapsize = 30*self.synthesized_beam()['bmaj'].to(u.mas).value
        fig.update_layout(coloraxis={'showscale': False, 'colorscale': 'Inferno'}, showlegend=False,
                          xaxis={'autorange': False, 'range': [mapsize, -mapsize]},
                          yaxis={'autorange': False, 'range': [-mapsize, mapsize]}, autosize=False)
        fig.update_xaxes(title_text="RA (mas)", constrain="domain")
        fig.update_yaxes(title_text="Dec (mas)", scaleanchor="x", scaleratio=1)
        # dirty_map_nat, laxis = obs.get_dirtymap(pixsize=1024, robust='natural', oversampling=4)
        # dirty_map_uni, laxis = obs.get_dirtymap(pixsize=1024, robust='uniform', oversampling=4)
        # fig1 = px.imshow(img=dirty_map_nat, x=laxis, y=laxis[::-1], labels={'x': 'RA (mas)', 'y': 'Dec (mas)'}, \
        #         aspect='equal')
        # fig2 = px.imshow(img=dirty_map_uni, x=laxis, y=laxis[::-1], labels={'x': 'RA (mas)', 'y': 'Dec (mas)'}, \
        #         aspect='equal')
        # fig = make_subplots(rows=1, cols=2, subplot_titles=('Natural weighting', 'Uniform weighting'),
        #                     shared_xaxes=True, shared_yaxes=True)
        # fig.add_trace(fig1.data[0], row=1, col=1)
        # fig.add_trace(fig2.data[0], row=1, col=2)
        # mapsize = 30*obs.synthesized_beam()['bmaj'].to(u.mas).value
        # fig.update_layout(coloraxis={'showscale': False, 'colorscale': 'Inferno'}, showlegend=False,
        #                   xaxis={'autorange': False, 'range': [mapsize, -mapsize]},
        #                   # This xaxis2 represents the xaxis for fig2.
        #                   xaxis2={'autorange': False, 'range': [mapsize, -mapsize]},
        #                   yaxis={'autorange': False, 'range': [-mapsize, mapsize]}, autosize=False)
        # fig.update_xaxes(title_text="RA (mas)", constrain="domain")
        # fig.update_yaxes(title_text="Dec (mas)", row=1, col=1, scaleanchor="x", scaleratio=1)

