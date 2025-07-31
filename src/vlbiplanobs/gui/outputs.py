# import io
import tempfile
from pathlib import Path
from typing import Optional, Literal
import numpy as np
from astropy import units as u
from dash import html, dcc
import dash_bootstrap_components as dbc
from borb import pdf
from vlbiplanobs import cli
from vlbiplanobs.gui import plots, inputs


def quantity2str(val: u.Quantity) -> str:
    return f"{val.value:.3g} {val.unit.to_string('unicode')}"


card = inputs.card


def card_result(number: str | list, label: str | list, id: str, extra_rows: Optional[list] = None,
                second_column_n: int = 0, second_column_content: Optional[list] = None) -> html.Div:
    return card(className='bg-primary opacity-10 p-0 m-0', style={'height': '100%', 'min-width': '200px'}, children=[
        html.Div(className='row', children=[
            html.Div(className=f'col-{12-second_column_n} text-start text-wrap: pretty', children=[
                html.Div(className='numbers', children=[
                    html.P(className='text-white text-sm mb-0 opacity-7 font-weight-bold',
                           children=label),
                    html.H3(className='text-white font-weight-bolder mb-0 justify-content-center',
                            id=f'{id}-value', children=number)
                ])
            ]),
            html.Div(className=f'col-{second_column_n} text-end',
                     children=second_column_content,
                     style={'align-content': 'center', 'text-wrap': 'pretty'})
            if second_column_n > 0 else '',
            html.Div(className='col-12', children=extra_rows) if extra_rows is not None else ''])])


def message_card(title: str | list, body: str | list, mode: Literal['danger', 'info', 'warning'],
                 icon: Optional[str] = None) -> html.Div:
    """Shows a card with a title and a body, with an icon and a mode that defines the color of the card:
    'info', 'warning', 'danger'. If 'icon' is not provided, it will take an icon by default depending on
    the mode.
    """
    assert mode in ('danger', 'info', 'warning')
    if icon is None:
        match mode:
            case 'danger':
                icon = 'fa fa-solid fa-triangle-exclamation'
            case 'info':
                icon = 'fa fa-solid fa-info'
            case 'warning':
                icon = 'fa fa-solid fa-exclamation'

    colors = {'danger': '#ef4444', 'warning': '#f4e2a8 !important', 'info': '#A6D4E8'}

    return html.Div(className='col-12 m-0 p-2', children=html.Div(
        className=f'card m-0 bg-{mode} opacity-10 border-radius-lg',
        children=html.Div(className='card-body p-2 position-relative', children=[
            html.Div(className='row align-items-center', children=[
                html.Div(className='col-2 pl-3 text-center',
                         children=html.Div(
                            html.I(className=f"{icon} text-xl text-{mode}-light",
                                   style={'font-size': '3rem',
                                          'color': colors[mode]}))),
                html.Div(className='col-10 text-start align-items-center', children=[
                    html.H5(className='text-white font-weight-bolder mb-0 mt-0',
                            children=title),
                    html.Span(className='text-white text-sm', children=body)])])])))


def warning_card(title: str | list, body: str | list, icon: Optional[str] = None) -> html.Div:
    return message_card(title, body, 'warning', icon)


def error_card(title: str | list, body: str | list, icon: Optional[str] = None) -> html.Div:
    return message_card(title, body, 'danger', icon)


def info_card(title: str | list, body: str | list, icon: Optional[str] = None) -> html.Div:
    return message_card(title, body, 'info', icon)


def warning_low_high_freq(o: Optional[cli.VLBIObs] = None) -> html.Div:
    """For observations at high frequencies, it shows a warning mentioning that
    phase referencing is not possible.
    For observations at low frequencies, it shows a warning mentioning that there
    may be significant RFI and the thermal rms may be affected.
    """
    if o is None:
        return html.Div()

    if o.frequency < 3*u.GHz:
        return warning_card("Significant RFI may be expected",
                            "At these low frequencies, a significant (~10%) of the data may be lost.",
                            icon='fa fa-solid fa-repeat')
    elif o.frequency > 70*u.GHz:
        return warning_card("Phase referencing is not possible",
                            "At these high frequencies, your target must be bright enough "
                            "to directly fringe on it.", icon='fa fa-solid fa-repeat')

    return html.Div()


def sun_warning(o: Optional[cli.VLBIObs] = None) -> html.Div:
    """Returns a warning card if the Sun is too close to the target source.
    Otherwise it returns an empty Div.
    """
    if o is None:
        return html.Div()

    sun_const = o.sun_constraint()
    sun_limit = o.sun_limiting_epochs()
    if not o.fixed_time:
        # if len(sun_limits := list(sun_limit.values())[0]) > 0:
        if list(sun_const.values())[0] is not None:
            sun_limits = list(sun_limit.values())[0]
            sun_const_src = sun_const[list(sun_const.keys())[0]]
            assert sun_const_src is not None, \
                "And error occured while checking if the Sun gets too close to the source."
            assert len(sun_limits) > 0, \
                "And error occured while checking if the Sun gets too close to the source."
            t0, t1 = sun_limits[0].datetime, sun_limits[-1].datetime
            if t0 == t1:
                return warning_card("The Sun gets too close to the source",
                                    f"On {t0.strftime('%-d %B')} "
                                    f"(with a minimum separation of {sun_const_src.value:.0f}"
                                    f"{sun_const_src.unit.to_string('unicode')}).", icon='fa fa-solid fa-sun')
            elif t0.month == t1.month:
                return warning_card("The Sun gets too close to the source",
                                    f"During {t0.day}-{t1.day} {t0.strftime('%B')} "
                                    f"(with a minimum separation of {sun_const_src.value:.0f}"
                                    f"{sun_const_src.unit.to_string('unicode')}).", icon='fa fa-solid fa-sun')
            elif t0.year == t1.year:
                return warning_card("The Sun gets too close to the source",
                                    f"From {t0.strftime('%-d %B')} to {t1.strftime('%-d %B')} "
                                    f"(with a minimum separation of {sun_const_src.value:.0f}"
                                    f"{sun_const_src.unit.to_string('unicode')}).", icon='fa fa-solid fa-sun')
            else:
                return warning_card("The Sun gets too close to the source",
                                    f"From {t0.strftime('%-d %B %Y')} to {t1.strftime('%-d %B %Y')} "
                                    f"(with a minimum separation of {sun_const_src.value:.0f}"
                                    f"{sun_const_src.unit.to_string('unicode')}).", icon='fa fa-solid fa-sun')
    else:
        if len(sun_separation := list(sun_const.values())) > 0 and sun_separation[0] is not None:
            return warning_card("The Sun is too close to the source!",
                                f"With a minimum separation of {sun_separation[0].value:.0f}"
                                f"{sun_separation[0].unit.to_string('unicode')} during the observation.",
                                icon='fa fa-solid fa-sun')

    return html.Div()


def ant_warning(o: Optional[cli.VLBIObs] = None) -> html.Div:
    """Returns a warning card if some antennas selected in the observation cannot observe.
    """
    if o is None:
        return html.Div()

    ants_excluded = []
    for src_is_visible in o.can_be_observed().values():
        ants_excluded = [ant for ant in src_is_visible if not src_is_visible[ant]]

    if not ants_excluded:
        return html.Div()

    return warning_card("Some antennas cannot observe the source",
                        inputs.parse_str_list([o.stations[ant].name for ant in ants_excluded]) +
                        ' will be out in the observation.', icon='fa fa-solid fa-satellite-dish')


def summary_freq_res(o: Optional[cli.VLBIObs] = None) -> html.Div:
    if o is None:
        return html.Div()

    try:
        chan_f = cli.optimal_units(o.bandwidth/(o.subbands*o.channels), [u.GHz, u.MHz, u.kHz, u.Hz])
        vel = cli.optimal_units((2.9979e5*u.km/u.s)*(o.bandwidth/(2*o.frequency)).decompose(),
                                [u.km/u.s, u.m/u.s])
        chan_v = cli.optimal_units(2*vel/(o.subbands*o.channels), [u.km/u.s])
    except Exception as e:
        # TODO: move it to a proper loggin
        print(f"Error computing frequency resolution: {e}")
        return html.Div()

    return card_result(f"{chan_f.value:.3g} {chan_f.unit.to_string('unicode')}",
                       "Freq. resolution", id='freq-res-frequency',
                       second_column_n=6, second_column_content=[
                            html.Div(className='numbers', children=[
                                html.P(className='text-white text-sm mb-0 opacity-7 font-weight-bold',
                                       children="Vel. resolution"),
                                html.H3(className='text-white font-weight-bolder mb-0 '
                                                  'justify-content-center',
                                        id='freq-res-velocity',
                                        children=f"{chan_v.value:.2g} "
                                                 f"{chan_v.unit.to_string('unicode')}")])
                       ], extra_rows=[
                            html.Br(),
                            html.Div(className='col-12', style={'color': 'var(--bs-gray-100)'},
                                     children=html.P(f"The total bandwidth of {o.bandwidth:n} is divided in "
                                                     f"{o.subbands} x {o.bandwidth/o.subbands:n} subbands, "
                                                     f"with {o.channels} spectral channels each.")
                                     if o.subbands > 1 else  # type: ignore
                                     html.P(f"The total bandwidth of {o.bandwidth:n} is recorded in "
                                            f"a single subband, with {o.channels} spectral channels."))])


def field_of_view(o: Optional[cli.VLBIObs] = None) -> html.Div:
    if o is None:
        return html.Div()

    if o.bandwidth_smearing() is None or o.time_smearing() is None:
        return html.Div()

    bw_smearing = cli.optimal_units(o.bandwidth_smearing(), [u.deg, u.arcmin, u.arcsec])
    tm_smearing = cli.optimal_units(o.time_smearing(), [u.deg, u.arcmin, u.arcsec])
    return card_result(f"{bw_smearing.value:.2n} {bw_smearing.unit.to_string('unicode')}",
                       "Bandwidth smearing", id='fov',
                       second_column_n=6, second_column_content=[
                            html.Div(className='numbers', children=[
                                html.P(className='text-white text-sm mb-0 opacity-7 font-weight-bold',
                                       children="Time smearing"),
                                html.H3(className='text-white font-weight-bolder mb-0 '
                                                  'justify-content-center',
                                        id=f'{id}-value',
                                        children=f"{tm_smearing.value:.2n} "
                                                 f"{tm_smearing.unit.to_string('unicode')}")])],
                       extra_rows=[
                            html.Br(),
                            html.Div(className='col-12', style={'color': 'var(--bs-gray-100)'},
                                     children=html.P("The field of view (FoV) will be limited to "
                                                     "this radius after correlation due to time "
                                                     "and frequency smearing, "
                                                     "considering a 10% loss."))])


def data_size(o: Optional[cli.VLBIObs] = None) -> html.Div:
    if o is None:
        return html.Div()

    ds = cli.optimal_units(o.datasize(), [u.TB, u.GB, u.MB])
    return card_result(f"{ds.value:.2n} {ds.unit.to_string('unicode')}",
                       "Data size", id='data-size',
                       extra_rows=[
                        html.Br(),
                        html.Div(className='col-12', style={'color': 'var(--bs-gray-100)'},
                                 children=html.P("This is the expected size that the correlated data will"
                                                 " have in FITS format. Note that during data reduction the required "
                                                 "space will be larger (e.g. about x3 for MS files)."))
                       ])


def obs_time(o: Optional[cli.VLBIObs] = None) -> html.Div:
    if o is None or o.ontarget_time is None:
        return html.Div()

    tt = cli.optimal_units(list(o.ontarget_time.values())[0], [u.h, u.min])
    dur = cli.optimal_units(o.duration, [u.h, u.min])
    wavelength = cli.optimal_units(o.wavelength, [u.m, u.cm, u.mm])
    return card_result(f"{dur.value:.3g} {dur.unit.to_string('unicode')}",
                       "Observing Time", id='obs-time-dur',
                       second_column_n=6, second_column_content=[
                            html.Div(className='numbers', children=[
                                html.P(className='text-white text-sm mb-0 opacity-7 font-weight-bold',
                                       children="Time on Source"),
                                html.H3(className='text-white font-weight-bolder mb-0 '
                                                  'justify-content-center',
                                        id='obs-time-target',
                                        children=f"{tt.value:.3g} "
                                                 f"{tt.unit.to_string('unicode')}")])
                       ], extra_rows=[
                            html.Br(),
                            html.Div(className='col-12', style={'color': 'var(--bs-gray-100)'},
                                     children=html.P("The planned observation will be conducted at a central "
                                        f"frequency of {o.frequency.to(u.GHz):.3n} "
                                        f"(wavelength of {wavelength:.3n})."))])

def rms(o: Optional[cli.VLBIObs] = None) -> html.Div:
    if o is None:
        return html.Div()

    try:
        if isinstance(thermal_noise := o.thermal_noise(), dict):
            rms = list(thermal_noise.values())[0]
        else:
            rms = thermal_noise

        out_rms: list[str] = [
            quantity2str(cli.optimal_units(rms, [u.MJy/u.beam, u.kJy/u.beam, u.Jy/u.beam,
                                                 u.mJy/u.beam, u.uJy/u.beam])),
            quantity2str(cli.optimal_units(rms/np.sqrt(1*u.min / (o.duration if o.duration is not None
                                                                  else 24*u.h)),
                                           [u.MJy/u.beam, u.kJy/u.beam, u.Jy/u.beam,
                                            u.mJy/u.beam, u.uJy/u.beam])),
            quantity2str(cli.optimal_units(rms*np.sqrt(o.subbands * o.channels),
                                           [u.MJy/u.beam, u.kJy/u.beam, u.Jy/u.beam,
                                            u.mJy/u.beam, u.uJy/u.beam]))]
    except Exception as e:
        # TODO: change this to a proper logging
        print(f"Error computing rms card: {e}")
        return html.Div()

    src = list(o.ontarget_time.keys())[0]
    return card_result(out_rms[0], f'rms thermal noise (for {o.ontarget_time[src]:.2g} on-target)',
                       id='rms', extra_rows=[html.Br(), html.Div(className='row', children=[
                            html.Div(className='col-6 text-start px-0 pb-2', children=[
                                html.H5(className='mb-0 font-weight-bolder text-light',
                                        id='rms-per-channel-value', children=out_rms[2]),
                                html.Label(className='text-sm mb-0 text-light',
                                           children="per spectral channel", id='tooltip-rms-channel')],
                                     style={'text-wrap': 'pretty'}),
                            dbc.Tooltip("Theoretical rms noise to obtain under ideal conditions "
                                        "in a single spectral channel when integrating over the whole "
                                        " observation.", target='tooltip-rms-channel'),
                            html.Div(className='col-6 text-end px-0', children=[
                                html.H5(className='mb-0 font-weight-bolder text-light',
                                        id='rms-per-time-value', children=out_rms[1]),
                                html.Label(className='text-sm mb-0 text-light',
                                           children="on 1-min integration", id='tooltip-rms-time'),
                                dbc.Tooltip("Theoretical rms noise to obtain under ideal conditions "
                                            "when integrating the data from the full bandwidth for "
                                            "one minute time integration.", target='tooltip-rms-time')],
                                     style={'text-wrap': 'pretty'})]),
                                   html.Div(className='row', children=[
                                       html.Div(html.Button('View sensitivity per baseline',
                                                            id='button-sensitivity-baseline',
                                                            className='btn btn-white btn-sm w-100 mb-0 active',
                                                            style={'position': 'sticky', 'top': '20px',
                                                                   'color': '#9DB7C4',
                                                                   'box-shadow': 'none'}))])])


def baseline_sensitivities(o: Optional[cli.VLBIObs] = None) -> html.Div:
    if o is None:
        return html.Div([dbc.ModalHeader(dbc.ModalTitle("Sensitivity per baseline")),
                         dbc.ModalBody("No information provided.")], id='sens-baseline-style')

    table_header = [html.Thead(html.Tr([html.Th("")] + [html.Th(s.codename,
                                                        className='text-left text-sm p-0 pb-2 ps-3')
                                                        for s in o.stations]))]
    all_baselines = np.unique(np.array([v.to(u.mJy/u.beam).value
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

        rows.append(html.Tr([html.Th(s.codename, className='text-left text-sm p-0 pb-2 ps-3'), *a_row]))

    table = dbc.Table(table_header + [html.Tbody(rows)], bordered=True,
                      className='table align-items-center mb-0')
    return html.Div([dbc.ModalHeader(
        dbc.ModalTitle("Sensitivity per baseline "
                       f"({(u.mJy/u.beam).to_string('unicode')}) for one-minute time integration")),
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
                       extra_rows=[html.Br(), html.Div(ellipse(bmaj="80px",
                                                               bmin=f"{bmin_elip}px",
                                                               pa=f"{-synth_beam['pa'].value+90:.1f}deg"),
                                                       className='row justift-content-center py-3',
                                                       style={'display': 'flex',
                                                              'justify-content': 'center'})])


def ellipse(bmaj, bmin, pa, color='white', z_index=1, position='relative', margin_top='', className=''):
    """Returns a html.Div element that draws an ellipse with a semimajor axis bmaj,
    semiminor axis bmin, and position angle (as defined in radio astronomy) pa.
    bmaj,bmin, pa must be strings recognized by HTML/CSS.
    """
    return html.Div(children=[],
                    style={'width': bmaj, 'height': bmin,
                           'border-radius': '50%', 'background': color, 'position': position,
                           'transform': f"rotate({pa})", 'z-index': z_index, 'justify-content': 'center',
                           'vertical-align': 'middle', 'margin-top': margin_top}, className=className)


def plot_elevations(o: Optional[cli.VLBIObs] = None) -> html.Div:
    return card([html.Div(className='card-header pb-0', children=html.H5('Source Elevation')),
                 dcc.Graph(id='fig-elevations', figure=plots.elevation_plot(o)),
                 dcc.Graph(id='fig-elevations2', figure=plots.elevation_plot_curves(o)),
                 html.Div(print_observability_ranges(o), id='out-elevations-info')])


def print_observability_ranges(o: Optional[cli.VLBIObs]) -> html.Div:
    if o is None:
        return html.Div()

    srcup = o.is_observable()

    ablockname = list(srcup.keys())[0]
    text: list = [html.Br()]
    when_everyone = o.when_is_observable(mandatory_stations='all',
                                         return_gst=not o.fixed_time)[ablockname]
    srcupalways = o.is_always_observable()
    if not when_everyone:
        text += ["The source cannot be observed by all stations at the same time."]
    else:
        if o.fixed_time:
            text += ["Everyone can observe the source on " +
                     ', '.join([t1.strftime('%d %b %Y %H:%M')+'-'+t2.strftime('%H:%M') +
                                ' UTC.' for t1, t2 in when_everyone])]
        else:
            text += ["Everyone can observe the source at " +
                     ', '.join([t1.to_string(sep=':', fields=2, pad=True) + '-' +
                                t2.to_string(sep=':', fields=2, pad=True) +
                                ' GST.' for t1, t2 in when_everyone])]

    text += [html.Br()]
    if any(srcupalways[ablockname].values()):
        ant_can = [ant for ant, b in srcupalways[ablockname].items() if b]
        if len(ant_can) <= len(o.stations)/2:
            text += [f"Only {', '.join(ant_can)} can observe the source at all times.", html.Br()]
        elif len(ant_can) < len(o.stations):
            text += ["All antennas can observe the source at all times except " +
                     ', '.join([ant for ant in o.stations.station_codenames if ant not in ant_can]),
                     '.', html.Br()]
    else:
        text += [f"{'And' if not when_everyone else 'But'} there are no antennas that can observe "
                 "the source at all time.", html.Br()]

    min_stat = 3 if len(o.stations) > 3 else min(2, len(o.stations))
    if len(o.stations) > 2:
        text += [f"The optimal visibility range (> {min_stat} antennas) occurs on "]
        if not o.fixed_time:
            text += [', '.join([t1.to_string(sep=':', fields=2, pad=True) + '-' +
                                (t2 + (24*u.hourangle if np.abs(t1 - t2) < 0.1*u.hourangle
                                       else 0.0*u.hourangle)).to_string(sep=':', fields=2, pad=True) +
                                ' GST.'
                                for t1, t2 in o.when_is_observable(min_stations=min_stat,
                                                                   return_gst=True)[ablockname]])]
        else:
            text += [', '.join([t1.strftime('%d %b %Y %H:%M')+'-'+t2.strftime('%H:%M') + ' UTC.'
                                for t1, t2 in o.when_is_observable(min_stations=min_stat)[ablockname]])]

    return html.Div(text)


def plot_uv_coverage(o: Optional[cli.VLBIObs] = None) -> html.Div:
    return card([html.Div(className='card-header pb-0', children=html.H5('(u, v) Coverage')),
                 dcc.Graph(id='fig-uv-coverage', figure=plots.uvplot(o)), html.P(""),
                 html.Label("Highlight antennas:", id='select-ant-uv-label',
                            htmlFor='select-ant-uv-plot'),
                 dcc.Dropdown(multi=True, id="select-antenna-uv-plot",
                              options=put_antenna_options(o)),  # type: ignore
                 html.Br(), html.Br(),
                 html.Div(className='row', id='out-uv-coverage-info', children=print_baseline_lengths(o))
                 ])


def put_antenna_options(o: Optional[cli.VLBIObs] = None) -> list:
    if o is None:
        return []

    return [{'label': f"{ant.name} ({ant.codename})", 'value': ant.codename} for ant in o.stations]


def print_baseline_lengths(o: Optional[cli.VLBIObs] = None) -> html.Div:
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

    return html.Div(className='row',
                    children=[html.Div(className='col-6 text-start text-wrap: pretty', children=[
                        html.Label(className='text-sm mb-0', children="Longest baseline"),
                        html.H5(className='mb-0 font-weight-bolder', id='baseline-long-value',
                                children=f"{cli.optimal_units(longest_bl[1], [u.km, u.m]):.5n} "
                                         f"({longest_bl_lambda.value:.3n} "
                                         f"{longest_bl_lambda.unit.name[0]}λ)"),
                        html.Label(className='text-sm mb-0',
                                   children=f"{o.stations[ant_l1].name} - "
                                            f"{o.stations[ant_l2].name} ({ant_l1}-{ant_l2})")]),
                              html.Div(className='col-6 text-end text-wrap: pretty', children=[
                                html.Label(className='text-sm mb-0', children="Shortest baseline"),
                                html.H5(className='mb-0 font-weight-bolder', id='baseline-short-value',
                                        children=f"{cli.optimal_units(shortest_bl[1], [u.km, u.m]):.5n} "
                                                 f"({shortest_bl_lambda.value:.3n} "
                                                 f"{shortest_bl_lambda.unit.name[0]}λ)"),
                                html.Label(className='text-sm mb-0',
                                           children=f"{o.stations[ant_s1].name} - "
                                                    f"{o.stations[ant_s2].name} ({ant_s1}-{ant_s2})")])])


def plot_worldmap(o: Optional[cli.VLBIObs] = None) -> html.Div:
    return card(html.Div(className='justify-content-center',
                         children=[dcc.Graph(figure=plots.plot_worldmap_stations(o),
                                             id='fig-worldmap',
                                             config={'showLink': False,
                                                     'displaylogo': False}),
                                   html.P("")]))


def download_button_div() -> html.Div:
    """Returns the placeholder for the button that allows to download the PDF with the card_results
    of the observation.
    """
    # return html.Div(html.Div(html.Div(), hidden=True, id='download-summary-div'),
    return html.Div(html.Div(download_button(), hidden=True, id='download-summary-div'),
                    className='col-6', style={'position': 'relative'})


def download_button() -> html.Div:
    """Returns the button to compute the observation
    """
    return html.Div([dbc.Spinner(id='downloading', color='#004990',
                                 children=html.Div(id='downloading-div')),
                     dbc.Button('Export Summary as PDF',
                                id='button-download', color='secondary', outline=True,
                                className='btn btn-lg btn-outline-secondary text-bolder '
                                            'mx-auto w-75 m-4 p-2',
                                style={'position': 'sticky', 'top': '20px'})],
                          className='d-flex align-items-center justify-content-center',
                          style={'gap': '5px'})


def summary_pdf(o: cli.VLBIObs, show_figure: bool = True):
    """Creates a PDF file with the summary of the observation and includes the elevation plot figure.
    """
    if o is None:
        raise ValueError("Observation cannot be None")

    doc: pdf.Document = pdf.Document()
    page = pdf.Page()
    doc.append_page(page)
    layout: pdf.PageLayout = pdf.SingleColumnLayout(page)
    layout.append_layout_element(pdf.Paragraph("EVN Observation Planner - Summary Report", font_size=20,
                             font='Helvetica-bold'))
                             # font='Helvetica-bold', horizontal_alignment=pdf.Alignment.CENTERED))
    text = f"Observation to be conducted at {o.band.replace('cm', ' cm')}"
    if o.fixed_time:
        if o.times[0].datetime.date() == o.times[-1].datetime.date():
            layout.append_layout_element(pdf.Paragraph(f"{text} from {o.times[0].strftime('%d %b %Y %H:%M')}–"
                                     f"{o.times[-1].strftime('%H:%M')} UTC."))
        elif (o.times[-1] - o.times[0]) < 24*u.h:
            layout.append_layout_element(pdf.Paragraph(f"{text} from {o.times[0].strftime('%d %b %Y %H:%M')}–"
                                     f"{o.times[-1].strftime('%H:%M')} (+1d) UTC."))
        else:
            layout.append_layout_element(pdf.Paragraph(f"{text} from {o.times[0].strftime('%d %b %Y %H:%M')} to "
                                     f"{o.times[-1].strftime('%d %b %H:%M')} UTC."))
    else:
        min_stat = 3 if len(o.stations) > 3 else min(2, len(o.stations))
        if len(o.stations) > 2:
            srcup = o.is_observable()
            if not srcup.items():
                layout.append_layout_element(pdf.Paragraph(f"{text}, but no epoch specified."))

            for ablockname, antbool in srcup.items():
                gst_range = (', '.join([t1.to_string(sep=':', fields=2, pad=True) + '--' +
                                        (t2 + (24*u.hourangle if np.abs(t1 - t2) < 0.1*u.hourangle
                                         else 0.0*u.hourangle)).to_string(sep=':', fields=2, pad=True) +
                                        ' GST.' for t1, t2
                                        in o.when_is_observable(min_stations=min_stat,
                                                                return_gst=not o.fixed_time)[ablockname]]))
                layout.append_layout_element(pdf.Paragraph(f"{text}. Optimal visibility range (> {min_stat} antennas) at "
                                         f"{gst_range}{' for '+ablockname if len(srcup) > 1 else ''}"))

    sun_const = o.sun_constraint()
    sun_limit = o.sun_limiting_epochs()
    for ablockname in o.sun_limiting_epochs():
        if not o.fixed_time:
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
                layout.append_layout_element(pdf.Paragraph(text, font_color=pdf.HexColor("#FF0000")))
        else:
            if (sun_const[ablockname] is not None) and not (not sun_const[ablockname]):
                text = "Note the the Sun is too close to the source"
                text += f" {ablockname}" if len(o.sun_limiting_epochs()) > 1 else " "
                text += f"({sun_const[ablockname]:.02f} away)."
                layout.append_layout_element(pdf.Paragraph(text, font_color=pdf.HexColor("#FF0000")))

    if o.duration is not None:
        layout.append_layout_element(pdf.Paragraph("With a total duration of "
                                 f"{cli.optimal_units(o.duration, [u.h, u.min, u.s]):.01f} "
                                 f"({cli.optimal_units(o.ontarget_time[list(o.ontarget_time.keys())[0]],
                                                       [u.h, u.min, u.s]):.01f} on target). "
                                 f"Total output FITS file size: {o.datasize():.2f}."))

    layout.append_layout_element(pdf.Paragraph(f"Participating stations ({len(o.stations)}): "
                             f"{', '.join(o.stations.station_codenames)}."))
    if not o.scans:
        layout.append_layout_element(pdf.Paragraph("No sources defined."))
    else:
        for ablock in o.scans.values():
            temp = '\n'.join([f"{s.name} ({s.coord.to_string('hmsdms')})." for s in ablock.sources()])
            layout.append_layout_element(pdf.Paragraph(f"Target source: {temp}"))

    # NOTE: for my future self: I do not like how units are displayed now, but using the unicode output
    # Makes the PDF writting to break because apparetly Helvetiva (and the other default fonts) do not
    # support symbols like "^-".  I tried!
    if None not in (o.datarate, o.bandwidth, o.subbands):
        val = cli.optimal_units(o.datarate, [u.Gbit/u.s, u.Mbit/u.s])
        layout.append_layout_element(pdf.Paragraph(f"\nData rate of {val:.0f}, "
                                 "producing a total bandwidth of "
                                 f"{cli.optimal_units(o.bandwidth, [u.MHz, u.GHz])}, "
                                 f" divided in {o.subbands} x {int(o.bandwidth.value/o.subbands)}-"
                                 f"{o.bandwidth.unit} subbands, with {o.channels} channels each, "
                                 f"{o.polarizations} polarization, and {o.inttime:.01f} integration time."))
    else:
        layout.append_layout_element(pdf.Paragraph("No setup (data rate, bandwidth, number of subbands) specified."))

    if not o.scans.values():
        rms = cli.optimal_units(o.thermal_noise(),
                                [u.Jy/u.beam, u.mJy/u.beam, u.uJy/u.beam])
        rms_chan = cli.optimal_units(rms/np.sqrt(1*u.min /
                                     (o.duration if o.duration is not None else 24*u.h)),
                                     [u.MJy/u.beam, u.kJy/u.beam, u.Jy/u.beam,
                                      u.mJy/u.beam, u.uJy/u.beam])
        rms_min = cli.optimal_units(rms*np.sqrt(o.subbands*o.channels),
                                    [u.MJy/u.beam, u.kJy/u.beam, u.Jy/u.beam,
                                     u.mJy/u.beam, u.uJy/u.beam])
        layout.append_layout_element(pdf.Paragraph(f"Thermal rms noise: "
                                 f"{rms:.3g}\n"
                                 f" ({rms_chan:.3g} per spectral "
                                 f"channel and {rms_min:.3g}"
                                 " per one-minute time integration."))
        if not o.sourcenames:
            synth_beam = o.synthesized_beam()[list(o.synthesized_beam().keys())[0]]
        else:
            synth_beam = o.synthesized_beam()[o.sourcenames[0]]

        bmaj = cli.optimal_units(synth_beam['bmaj'], [u.deg, u.arcmin, u.arcsec, u.mas, u.uas])
        bmin = synth_beam['bmin'].to(bmaj.unit)
        layout.append_layout_element(pdf.Paragraph(f"Synthesized beam (approx for a random source): "
                                 f"{bmaj.value:2.1f} x {bmin:2.1f}"
                                 f", {synth_beam['pa'].value:2.0f}º."))

    else:
        for ablock in o.scans.values():
            for src in o.sources():
                if len(o.scans) > 1:
                    layout.append_layout_element(pdf.Paragraph(f"For the source {src.name}", font='Helvetica-bold'))

                rms = cli.optimal_units(o.thermal_noise()[src.name],
                                        [u.Jy/u.beam, u.mJy/u.beam, u.uJy/u.beam])
                rms_chan = cli.optimal_units(rms/np.sqrt(1*u.min /
                                             (o.duration if o.duration is not None else 24*u.h)),
                                             [u.MJy/u.beam, u.kJy/u.beam, u.Jy/u.beam,
                                              u.mJy/u.beam, u.uJy/u.beam])
                rms_min = cli.optimal_units(rms*np.sqrt(o.subbands*o.channels),
                                            [u.MJy/u.beam, u.kJy/u.beam, u.Jy/u.beam,
                                             u.mJy/u.beam, u.uJy/u.beam])
                layout.append_layout_element(pdf.Paragraph(f"Thermal rms noise: "
                                         f"{rms:.3g}\n"
                                         f" ({rms_chan:.3g} per spectral "
                                         f"channel and {rms_min:.3g}"
                                         " per one-minute time integration."))
                synth_beam = o.synthesized_beam()[src.name]
                bmaj = cli.optimal_units(synth_beam['bmaj'], [u.deg, u.arcmin, u.arcsec, u.mas, u.uas])
                bmin = synth_beam['bmin'].to(bmaj.unit)
                temp = f" for {src.name}" if len(o.scans.values()) > 1 else ""
                layout.append_layout_element(pdf.Paragraph(f"Synthesized beam{temp}: "
                                         f"{bmaj.value:2.1f} x {bmin:2.1f}"
                                         f", {synth_beam['pa'].value:2.0f}º."))

    bw_smearing = cli.optimal_units(o.bandwidth_smearing(), [u.deg, u.arcmin, u.arcsec])
    tm_smearing = cli.optimal_units(o.time_smearing(), [u.deg, u.arcmin, u.arcsec])
    layout.append_layout_element(pdf.Paragraph(f"Field of view limited to {bw_smearing:.2g} (from frequency smearing) "
                             f"and {tm_smearing:.2g} (from time smearing), considering 10% loss."))

    if len(o.scans) > 0 and show_figure:
        fig = plots.elevation_plot(o, show_colorbar=True)
        assert fig is not None, "An image could not be created for the PDF"
        with tempfile.NamedTemporaryFile(suffix='.jpg', delete=False) as tempfig:
            fig.update_layout(plot_bgcolor='white', paper_bgcolor='white', font_color='black')
            fig.write_image(tempfig.name, scale=2, width=800)
            figpath = Path(tempfig.name)
            print(tempfig.name)

        # layout.append_layout_element(pdf.Image(figpath, width=414, height=265, horizontal_alignment=pdf.Alignment.CENTERED))
        layout.append_layout_element(pdf.Image(figpath, size=(414, 265)))

    tmp = tempfile.NamedTemporaryFile(suffix='.pdf', delete=False)
    pdf.PDF.write(where_to=tmp.name, what=doc)
    print(f"File at {tmp.name}.")
        # buffer.seek(0)

    return tmp.name

    # def get_fig_dirty_map(self):
    #     raise NotImplementedError
    #     # Right now I only leave the natural weighting map (the uniform does not always correspond to the true one)
    #     dirty_map_nat, laxis = self.get_dirtymap(pixsize=1024, robust='natural', oversampling=4)
    #     # TODO: uncomment these two lines and move them outside observation. Flexibility
    #     # fig1 = px.imshow(img=dirty_map_nat, x=laxis, y=laxis[::-1], labels={'x': 'RA (mas)', 'y': 'Dec (mas)'}, \
    #     #         aspect='equal')
    #     # fig = make_subplots(rows=1, cols=1, subplot_titles=('Natural weighting',), shared_xaxes=True, shared_yaxes=True)
    #     fig.add_trace(fig1.data[0], row=1, col=1)
    #     mapsize = 30*self.synthesized_beam()['bmaj'].to(u.mas).value
    #     fig.update_layout(coloraxis={'showscale': False, 'colorscale': 'Inferno'}, showlegend=False,
    #                       xaxis={'autorange': False, 'range': [mapsize, -mapsize]},
    #                       yaxis={'autorange': False, 'range': [-mapsize, mapsize]}, autosize=False)
    #     fig.update_xaxes(title_text="RA (mas)", constrain="domain")
    #     fig.update_yaxes(title_text="Dec (mas)", scaleanchor="x", scaleratio=1)
    #     # dirty_map_nat, laxis = obs.get_dirtymap(pixsize=1024, robust='natural', oversampling=4)
    #     # dirty_map_uni, laxis = obs.get_dirtymap(pixsize=1024, robust='uniform', oversampling=4)
    #     # fig1 = px.imshow(img=dirty_map_nat, x=laxis, y=laxis[::-1], labels={'x': 'RA (mas)', 'y': 'Dec (mas)'}, \
    #     #         aspect='equal')
    #     # fig2 = px.imshow(img=dirty_map_uni, x=laxis, y=laxis[::-1], labels={'x': 'RA (mas)', 'y': 'Dec (mas)'}, \
    #     #         aspect='equal')
    #     # fig = make_subplots(rows=1, cols=2, subplot_titles=('Natural weighting', 'Uniform weighting'),
    #     #                     shared_xaxes=True, shared_yaxes=True)
    #     # fig.add_trace(fig1.data[0], row=1, col=1)
    #     # fig.add_trace(fig2.data[0], row=1, col=2)
    #     # mapsize = 30*obs.synthesized_beam()['bmaj'].to(u.mas).value
    #     # fig.update_layout(coloraxis={'showscale': False, 'colorscale': 'Inferno'}, showlegend=False,
    #     #                   xaxis={'autorange': False, 'range': [mapsize, -mapsize]},
    #     #                   # This xaxis2 represents the xaxis for fig2.
    #     #                   xaxis2={'autorange': False, 'range': [mapsize, -mapsize]},
    #     #                   yaxis={'autorange': False, 'range': [-mapsize, mapsize]}, autosize=False)
    #     # fig.update_xaxes(title_text="RA (mas)", constrain="domain")
    #     # fig.update_yaxes(title_text="Dec (mas)", row=1, col=1, scaleanchor="x", scaleratio=1)

