from collections import defaultdict
from typing import Optional
from datetime import datetime
import math
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.colors import get_colorscale
from vlbiplanobs import sources


_MINOR_TICKS = dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)
_TRANSPARENT = 'rgba(0,0,0,0)'
_VIRIDIS = get_colorscale('Viridis')


def _figure(traces: list[dict], layout: dict) -> go.Figure:
    """Build a figure from plain trace/layout dicts, skipping plotly's per-property validation.

    Creating graph objects and calling update_layout validates (and parses the name of) every single
    property, which takes several times longer than computing the data that is plotted. The figures of
    the web app are rebuilt on every change of the inputs, so they are assembled as plain dicts instead.
    The dicts must then be written in plotly's canonical form: no "magic underscore" keys (xaxis_title)
    and no shorthands (title='...' must be title=dict(text='...')).

    Parameters
    ----------
    traces : list[dict]
        Traces of the figure, each one including its 'type'.
    layout : dict
        Layout of the figure.

    Returns
    -------
    go.Figure
        Plotly figure.
    """
    try:
        return go.Figure(data=traces, layout=layout, _validate=False)
    except (TypeError, ValueError):  # plotly version without the (semi-private) _validate option
        return go.Figure(data=traces, layout=layout)


def _compute_gst_ticks(times_dt: list, gst_hours: list) -> tuple[list, list]:
    """Compute rounded major-tick positions for a GST top axis paired with a UTC bottom axis.

    The GST values are unwrapped to handle midnight wraparound, then rounded tick
    positions are mapped back to UTC datetimes via linear interpolation.

    Parameters
    ----------
    times_dt : list
        Sequence of UTC datetime objects (length n, monotonic).
    gst_hours : list
        Sequence of GST hours in [0, 24) (length n, same sampling as times_dt).

    Returns
    -------
    tuple[list, list]
        (tickvals, ticktext) where tickvals is a list of UTC datetime objects
        placed on rounded GST values and ticktext is the matching list of "HH:MM" strings.
    """
    n = len(times_dt)
    if n < 2 or len(gst_hours) != n:
        return [], []

    unwrapped = [float(gst_hours[0])]
    for i in range(1, n):
        d = float(gst_hours[i]) - float(gst_hours[i - 1])
        if d < -12.0:
            d += 24.0
        elif d > 12.0:
            d -= 24.0
        unwrapped.append(unwrapped[-1] + d)

    g_start, g_end = unwrapped[0], unwrapped[-1]
    span = g_end - g_start
    if span <= 0:
        return [], []

    if span > 18.0:
        step = 3.0
    elif span > 12.0:
        step = 2.0
    elif span > 6.0:
        step = 1.0
    elif span > 3.0:
        step = 0.5
    elif span > 1.5:
        step = 0.25
    elif span > 0.5:
        step = 1.0 / 6.0
    else:
        step = 1.0 / 12.0

    unwrapped_arr = np.asarray(unwrapped)
    eps = step * 1e-6
    first = math.ceil((g_start - eps) / step) * step

    tickvals: list = []
    ticktext: list = []
    n_steps = int(math.floor((g_end - first) / step + eps)) + 1
    for k in range(max(n_steps, 0)):
        g = first + k * step
        if g < g_start - eps or g > g_end + eps:
            continue
        idx = int(np.searchsorted(unwrapped_arr, g) - 1)
        if idx < 0:
            idx = 0
        if idx >= n - 1:
            idx = n - 2
        denom = unwrapped_arr[idx + 1] - unwrapped_arr[idx]
        frac = 0.0 if denom == 0 else (g - unwrapped_arr[idx]) / denom
        t = times_dt[idx] + (times_dt[idx + 1] - times_dt[idx]) * frac

        gh = g % 24.0
        h = int(gh)
        m = int(round((gh - h) * 60.0))
        if m == 60:
            h += 1
            m = 0
        h = h % 24
        tickvals.append(t)
        ticktext.append(f"{h:02d}:{m:02d}")
    return tickvals, ticktext


def _build_gst_axis_config(times_dt: list, gst_hours: list, time_range: list) -> Optional[dict]:
    """Build the layout config for a GST top axis overlaying the UTC bottom axis.

    Returns None when the GST axis cannot be built (insufficient samples, no span).
    The returned dict contains the xaxis2 layout entry with rounded "HH:MM" ticks
    and an explicit range covering the same UTC window as the bottom axis.

    Parameters
    ----------
    times_dt : list
        Sequence of UTC datetime objects.
    gst_hours : list
        Sequence of GST hours in [0, 24).
    time_range : list
        UTC time range for the plot.

    Returns
    -------
    dict or None
        Layout config for the GST axis, or None if it cannot be built.

    Notes
    -----
    Does not use ``matches='x'`` / ``anchor='y'`` because they conflict with
    ``overlaying='x'`` in several plotly versions (the secondary axis silently fails
    to render).
    """
    tickvals, ticktext = _compute_gst_ticks(times_dt, gst_hours)
    if not tickvals:
        return None
    return dict(
        overlaying='x', side='top', type='date',
        range=time_range, autorange=False,
        tickmode='array', tickvals=tickvals, ticktext=ticktext,
        showline=True, linecolor='black', linewidth=1, mirror=False,
        ticks='inside', showgrid=False, zeroline=False,
        title=dict(text='Time (GST)'),
    )


def _apply_gst_top_axis(fig: go.Figure, axis_cfg: dict, y_anchor: float = 0.0) -> None:
    """Add the GST top axis (xaxis2) to fig plus a fully transparent anchor trace.

    Plotly only draws a secondary axis when at least one trace is bound to it and
    that trace contributes points. This helper plots two markers at y=y_anchor
    with fully transparent markers and disabled hover, so the axis is rendered
    without visually touching the data.

    Parameters
    ----------
    fig : go.Figure
        Plotly figure to add the axis to.
    axis_cfg : dict
        Axis configuration dictionary.
    y_anchor : float, optional
        Y-coordinate for the anchor trace. Default is 0.0.
    """
    cfg, anchor_trace = _gst_top_axis(axis_cfg, y_anchor)
    if cfg is None:
        return

    fig.update_layout(xaxis2=cfg)
    fig.add_trace(anchor_trace)


def _gst_top_axis(axis_cfg: dict, y_anchor: float = 0.0) -> tuple[Optional[dict], Optional[dict]]:
    """Return the GST top axis (xaxis2) layout and its transparent anchor trace, as plain dicts.

    Parameters
    ----------
    axis_cfg : dict
        Axis configuration dictionary.
    y_anchor : float, optional
        Y-coordinate for the anchor trace. Default is 0.0.

    Returns
    -------
    tuple[dict or None, dict or None]
        (xaxis2 layout, anchor trace), or (None, None) when there is no GST top axis.
    """
    cfg = axis_cfg.get('xaxis2_config')
    if cfg is None:
        return None, None

    rng = axis_cfg['xaxis_range']
    return cfg, dict(type='scatter', x=[rng[0], rng[-1]], y=[y_anchor, y_anchor], xaxis='x2', yaxis='y',
                     mode='markers', marker=dict(opacity=0, size=1), opacity=0,
                     showlegend=False, hoverinfo='skip')


def _get_axis_config(o):
    """Return axis configuration for elevation plots based on observation time mode.

    Parameters
    ----------
    o : VLBIObs
        VLBI observation object with times and gstimes attributes.

    Returns
    -------
    dict
        Dictionary with keys:
        - xaxis_range: [start, end] datetime range for the bottom (UTC) x-axis.
        - xaxis_title: title for the bottom x-axis.
        - xaxis2_config: layout dict for the secondary GST x-axis (None when not applicable).
    """
    times_dt = o.times.datetime
    time_range = [times_dt[0], times_dt[-1]]

    if o.fixed_time:
        gst_hours = [float(g.hour) for g in o.gstimes]
        xaxis2_config = _build_gst_axis_config(list(times_dt), gst_hours, time_range)
        return {
            'xaxis_range': time_range,
            'xaxis_title': 'Time (UTC)',
            'xaxis2_config': xaxis2_config,
        }
    return {'xaxis_range': time_range, 'xaxis_title': 'Time (GST)', 'xaxis2_config': None}


def elevation_plot_curves(o) -> Optional[go.Figure]:
    """Create the plot showing when antennas can observe a source (old style).

    Shows elevation on the y-axis.

    Parameters
    ----------
    o : VLBIObs
        VLBI observation object.

    Returns
    -------
    go.Figure or None
        Plotly figure, or None if observation is invalid.
    """
    if o is None or not o.scans:
        return None

    srcup = o.is_observable()
    elevs = o.elevations()
    times_dt = o.times.datetime
    traces: list[dict] = []
    for src_block in srcup:
        targets = o.scans[src_block].sources()  # sources.SourceType.TARGET)
        for ant in srcup[src_block]:
            y = np.full(len(times_dt), np.nan)
            y[srcup[src_block][ant]] = elevs[targets[0].name][ant][srcup[src_block][ant]].value
            traces.append(dict(type='scatter', x=times_dt, y=y, mode='lines', connectgaps=False, hoverinfo='none',
                               hovertemplate=f"<b>{o.stations[ant].name} ({o.stations[ant].codename})</b><br>"
                                             "<b>Time</b>: %{x}",  # .strftime('%H:%M')}",
                               name=o.stations[ant].name))

    # Get axis configuration based on observation time mode
    axis_cfg = _get_axis_config(o)
    xaxis2, anchor_trace = _gst_top_axis(axis_cfg)

    # When the bottom axis is UTC (fixed_time), the top axis mirrors it as GST,
    # so we keep the bottom axis in 'allticks' mirror mode only without a top axis.
    bottom_mirror = 'allticks' if xaxis2 is None else False

    # Gray bands marking the low elevations (below 20 and 10 degrees)
    low_elevation = dict(type='rect', xref='paper', yref='y', x0=0, y0=0, x1=1, fillcolor='gray', opacity=0.2,
                         layer='below', line=dict(width=0))
    layout = dict(shapes=[dict(low_elevation, y1=20), dict(low_elevation, y1=10)],
                  showlegend=True, hovermode='closest', title=dict(text=''), paper_bgcolor=_TRANSPARENT,
                  plot_bgcolor=_TRANSPARENT,
                  xaxis=dict(anchor='y', domain=[0.0, 1.0], title=dict(text=axis_cfg['xaxis_title']),
                             type="date", tickformat="%H:%M", showline=True, linecolor='black', linewidth=1,
                             mirror=bottom_mirror, ticks='inside', tickmode='auto', range=axis_cfg['xaxis_range'],
                             showgrid=False, zeroline=False, minor=_MINOR_TICKS),
                  yaxis=dict(anchor='x', domain=[0.0, 1.0], title=dict(text="Elevation (degrees)"),
                             showline=True, linecolor='black', linewidth=1, mirror='allticks',
                             ticks='inside', tickmode='auto', range=[0, 90],
                             showgrid=False, zeroline=False, minor=_MINOR_TICKS),
                  margin=dict(l=2, r=2, t=45, b=0),
                  legend=dict(x=0.01, y=0.99, xanchor='left', yanchor='top',
                              bgcolor='rgba(255, 255, 255, 0.7)',
                              bordercolor='rgba(0, 0, 0, 0.3)', borderwidth=1))
    if xaxis2 is not None:
        layout['xaxis2'] = xaxis2
        traces.append(anchor_trace)

    return _figure(traces, layout)


def elevation_plot(o, show_colorbar: bool = False) -> Optional[go.Figure]:
    """Create the plot showing when antennas can observe a source.

    Optimized version using Heatmap instead of individual traces.

    Parameters
    ----------
    o : VLBIObs
        VLBI observation object.
    show_colorbar : bool, optional
        Whether to show the colorbar. Default is False.

    Returns
    -------
    go.Figure or None
        Plotly figure, or None if observation is invalid.
    """
    if o is None or not o.scans:
        return None

    srcup = o.is_observable()
    elevs = o.elevations()
    times_dt = o.times.datetime
    time_labels = [t.strftime('%H:%M') for t in times_dt]
    n_times = len(times_dt)

    # One heatmap per block, with its antenna names so each subplot gets its own y-axis labels
    heatmaps: list[tuple[dict, list[str]]] = []
    for src_block in srcup:
        ant_names = list(srcup[src_block].keys())
        n_ants = len(ant_names)
        targets = o.scans[src_block].sources(sources.SourceType.TARGET)
        target_name = targets[0].name if targets else o.scans[src_block].sources()[0].name

        # Build elevation matrix (antennas x times)
        z_matrix = np.full((n_ants, n_times), np.nan)
        for anti, ant in enumerate(ant_names):
            visibility = np.array(srcup[src_block][ant], dtype=bool)
            z_matrix[n_ants - 1 - anti, visibility] = elevs[target_name][ant].value[visibility]

        # Create hover text matrix
        hover_text = []
        for ai in range(n_ants):
            station_name = o.stations[ant_names[n_ants - 1 - ai]].name
            hover_text.append([f"<b>{station_name}</b><br>Elevation: {z:.0f}º<br>Time: {label}"
                               if not np.isnan(z) else ""
                               for z, label in zip(z_matrix[ai].tolist(), time_labels)])

        heatmaps.append((dict(type='heatmap', x=times_dt, y=list(range(1, n_ants + 1)), z=z_matrix,
                              colorscale=_VIRIDIS, zmin=5, zmax=90, showscale=show_colorbar,
                              hoverinfo='text', text=hover_text, xgap=0, ygap=1), ant_names))

    # Get axis configuration based on observation time mode
    axis_cfg = _get_axis_config(o)
    # Only attach a GST top axis when the heatmap has a single subplot. With multiple
    # subplots, make_subplots already populates xaxis2/3/... so adding our overlay would
    # clash with their layout.
    single_subplot = len(srcup) == 1
    if not single_subplot:
        axis_cfg = {**axis_cfg, 'xaxis2_config': None}

    bottom_mirror = 'allticks' if axis_cfg.get('xaxis2_config') is None else False
    xaxis = dict(type="date", tickformat="%H:%M", showline=True, linecolor='black', linewidth=1,
                 mirror=bottom_mirror, ticks='inside', tickmode='auto', range=axis_cfg['xaxis_range'],
                 showgrid=False, zeroline=False, minor=_MINOR_TICKS)
    yaxis = dict(showline=True, linecolor='black', linewidth=1, mirror='allticks', ticks='inside',
                 showgrid=False, zeroline=False, minor=_MINOR_TICKS)
    layout = dict(showlegend=False, hovermode='closest', title=dict(text=''), paper_bgcolor=_TRANSPARENT,
                  plot_bgcolor=_TRANSPARENT, margin=dict(l=2, r=2, t=45, b=0))
    if show_colorbar:
        layout['coloraxis'] = dict(colorscale=_VIRIDIS, colorbar=dict(title=dict(text='Elevation (degrees)')))

    if single_subplot:
        # The usual case (always in the web app): a single plot, assembled without plotly's validation.
        heatmap, names = heatmaps[0]
        traces = [dict(heatmap, xaxis='x', yaxis='y')]
        layout['xaxis'] = dict(xaxis, anchor='y', domain=[0.0, 1.0], title=dict(text=axis_cfg['xaxis_title']))
        layout['yaxis'] = dict(yaxis, anchor='x', domain=[0.0, 1.0], title=dict(text="Antennas"),
                               tickmode='array', tickvals=list(range(1, len(names) + 1)),
                               ticktext=list(reversed(names)))
        xaxis2, anchor_trace = _gst_top_axis(axis_cfg)
        if xaxis2 is not None:
            layout['xaxis2'] = xaxis2
            traces.append(anchor_trace)

        return _figure(traces, layout)

    fig = make_subplots(rows=min([len(srcup), 4]), cols=len(srcup) // 4 + 1,
                        subplot_titles=[f"Elevations for {src_block}" for src_block in srcup])
    for src_i, (heatmap, _) in enumerate(heatmaps):
        fig.add_trace(go.Heatmap(**{key: value for key, value in heatmap.items() if key != 'type'}),
                      row=src_i % 4 + 1, col=src_i // 4 + 1)

    fig.update_layout(xaxis=dict(xaxis, title=dict(text=axis_cfg['xaxis_title'])),
                      yaxis=dict(title=dict(text="Antennas")), **layout)
    # Style all y axes, then label each subplot with its own block's antenna names.
    fig.update_yaxes(**yaxis)
    for src_i, (_, names) in enumerate(heatmaps):
        fig.update_yaxes(tickmode='array', tickvals=list(range(1, len(names) + 1)),
                         ticktext=list(reversed(names)), row=src_i % 4 + 1, col=src_i // 4 + 1)

    return fig


# uv-plot styling. The clientside highlight callback in callbacks.py receives these same constants (embedded as
# JSON) so that the Python-built figure and the browser-restyled figure always look identical.
UV_HIGHLIGHT_COLORS: list[str] = ["#FF0000", "#0000FF", "#008000", "#FFA500", "#800080",
                                  "#00FFFF", "#FF00FF", "#FFD700", "#00FF00", "#A52A2A",
                                  "#FFC0CB", "#808000", "#008080", "#000080", "#FF7F50",
                                  "#4B0082", "#FF8C00", "#40E0D0", "#6A5ACD", "#006400"]
UV_BASE_COLOR: str = 'black'
UV_BASE_SIZE: int = 2
UV_HIGHLIGHT_SIZE: int = 4


def baseline_has_antenna(baseline: str, antenna: str) -> bool:
    """Return True if antenna is one of the two stations of a baseline key.

    Baseline keys are built in observation.py as f"{codename1}-{codename2}" (station codenames
    never contain '-'). An exact match on the split parts is required so that e.g. 'Me' does not
    match 'Me1-Ef' and 'Wb' does not match 'Wb14-Ef'.

    Parameters
    ----------
    baseline : str
        Baseline key, e.g. 'Ef-Wb'.
    antenna : str
        Station codename to look for.

    Returns
    -------
    bool
        True if antenna is either end of the baseline.
    """
    return antenna in baseline.split('-')


def uv_marker_style(baseline: str, filter_antennas: Optional[list[str]] = None) -> dict:
    """Return the marker style of one baseline trace in the uv plot.

    The first antenna in filter_antennas that belongs to the baseline decides the highlight color
    (UV_HIGHLIGHT_COLORS cycled by that antenna's position). Must stay in sync with the JavaScript
    `uvStyleFor` function in callbacks.py.

    Parameters
    ----------
    baseline : str
        Baseline key, e.g. 'Ef-Wb'.
    filter_antennas : list[str] or None, optional
        Antennas to highlight. Default is None (no highlight).

    Returns
    -------
    dict
        {'color': str, 'size': int, 'highlighted': bool}.
    """
    for i, ant in enumerate(filter_antennas or []):
        if baseline_has_antenna(baseline, ant):
            return {'color': UV_HIGHLIGHT_COLORS[i % len(UV_HIGHLIGHT_COLORS)], 'size': UV_HIGHLIGHT_SIZE,
                    'highlighted': True}
    return {'color': UV_BASE_COLOR, 'size': UV_BASE_SIZE, 'highlighted': False}


def uvplot_from_baselines(bl_uv: dict, filter_antennas: Optional[list[str]] = None) -> go.Figure:
    """Create the uv-coverage figure with one Scattergl trace per baseline.

    Each trace is named after its baseline (e.g. 'Ef-Wb') and contains both the +uv and -uv points as
    float32 numpy arrays (plotly serializes them as compact base64 typed arrays). The hover label is the
    trace name via hovertemplate, so no per-point text is shipped. Highlighted traces are placed last so
    they are drawn on top.

    Parameters
    ----------
    bl_uv : dict
        Mapping baseline key -> array-like (or astropy Quantity) of shape (N, >=2) with u, v in lambda.
    filter_antennas : list[str] or None, optional
        Antennas to highlight in the plot. Default is None.

    Returns
    -------
    go.Figure
        Plotly figure (possibly with no traces if bl_uv is empty).
    """
    base_traces: list = []
    highlighted_traces: list = []
    for baseline, arr in bl_uv.items():
        values = np.asarray(arr.value if hasattr(arr, 'value') else arr)
        if values.ndim != 2 or values.shape[0] == 0:
            continue
        uu = values[:, 0].astype(np.float32)
        vv = values[:, 1].astype(np.float32)
        style = uv_marker_style(baseline, filter_antennas)
        trace = dict(type='scattergl', x=np.concatenate([uu, -uu]), y=np.concatenate([vv, -vv]), name=baseline,
                     mode='markers', marker=dict(color=style['color'], size=style['size']),
                     hovertemplate='%{fullData.name}<extra></extra>', showlegend=False)
        (highlighted_traces if style['highlighted'] else base_traces).append(trace)

    axis = dict(showline=True, linecolor='black', linewidth=1, mirror='allticks', ticks='inside',
                tickmode='auto', showgrid=False, zeroline=False, minor=_MINOR_TICKS)
    return _figure(base_traces + highlighted_traces, dict(
        showlegend=False,
        hovermode='closest',
        uirevision='uv',  # keep zoom/pan when the highlight callback restyles the figure
        title=dict(text=''),
        paper_bgcolor=_TRANSPARENT,
        plot_bgcolor=_TRANSPARENT,
        xaxis=dict(axis, title=dict(text='u (λ)'), constrain='domain'),
        yaxis=dict(axis, title=dict(text='v (λ)'), scaleanchor="x", scaleratio=1),
        margin=dict(l=0, r=0, t=0, b=0)))


def uvplot(o, filter_antennas: Optional[list[str]] = None) -> Optional[go.Figure]:
    """Create the uv-coverage figure (one trace per baseline) for the first source of an observation.

    Parameters
    ----------
    o : VLBIObs
        VLBI observation object.
    filter_antennas : list[str] or None, optional
        Antennas to highlight in the plot. Default is None.

    Returns
    -------
    go.Figure or None
        Plotly figure, or None if the observation has no scans or no uv data.
    """
    if o is None or not o.scans:
        return None

    bl_uv = o.get_uv_data()
    if not bl_uv:
        return None
    return uvplot_from_baselines(bl_uv[list(bl_uv.keys())[0]], filter_antennas)


def serialize_elevation_data(o) -> Optional[dict]:
    """Serialize elevation/observability data for deferred plot rendering.

    Returns a JSON-serializable dict with all data needed by elevation_plot_from_data
    and elevation_curves_from_data.

    Parameters
    ----------
    o : VLBIObs
        VLBI observation object.

    Returns
    -------
    dict or None
        Serialized elevation data, or None if no scans.
    """
    if o is None or not o.scans:
        return None

    srcup = o.is_observable()
    elevs = o.elevations()
    times_iso = [t.isoformat() for t in o.times.datetime]
    fixed_time = o.fixed_time
    gstimes_hours = [float(g.hour) for g in o.gstimes] if hasattr(o, 'gstimes') else []
    first_date_iso = o.times.datetime[0].date().isoformat()

    blocks = {}
    for src_block in srcup:
        ant_names = list(srcup[src_block].keys())
        targets = o.scans[src_block].sources(sources.SourceType.TARGET)
        target_name = targets[0].name if targets else o.scans[src_block].sources()[0].name

        observability = {ant: srcup[src_block][ant].astype(bool).tolist() for ant in ant_names}
        elevation_vals = {ant: elevs[target_name][ant].value.tolist() for ant in ant_names}

        blocks[src_block] = {
            'ant_names': ant_names,
            'observability': observability,
            'elevations': elevation_vals,
        }

    station_info = {ant.codename: {'name': ant.name, 'codename': ant.codename} for ant in o.stations}

    return {
        'times_iso': times_iso,
        'fixed_time': fixed_time,
        'gstimes_hours': gstimes_hours,
        'first_date_iso': first_date_iso,
        'blocks': blocks,
        'station_info': station_info,
    }


def _get_axis_config_from_data(data: dict) -> dict:
    """Build axis config from serialized elevation data.

    Mirrors _get_axis_config but reads from JSON-serialized payload produced by
    serialize_elevation_data. The GST top axis is only built for fixed-epoch
    observations that carry GST samples matching the time grid.

    Parameters
    ----------
    data : dict
        Serialized elevation data from serialize_elevation_data.

    Returns
    -------
    dict
        Axis configuration dictionary.
    """
    times = [datetime.fromisoformat(t) for t in data['times_iso']]
    time_range = [times[0], times[-1]]
    if data['fixed_time']:
        gst_hours = data.get('gstimes_hours') or []
        xaxis2_config = (_build_gst_axis_config(times, gst_hours, time_range)
                         if len(gst_hours) == len(times) else None)
        return {'xaxis_range': time_range, 'xaxis_title': 'Time (UTC)',
                'xaxis2_config': xaxis2_config}
    return {'xaxis_range': time_range, 'xaxis_title': 'Time (GST)', 'xaxis2_config': None}


def elevation_plot_from_data(data: dict, show_colorbar: bool = False) -> Optional[go.Figure]:
    """Render the heatmap elevation plot from serialized data.

    Parameters
    ----------
    data : dict
        Serialized elevation data from serialize_elevation_data.
    show_colorbar : bool, optional
        Whether to show the colorbar. Default is False.

    Returns
    -------
    go.Figure or None
        Plotly figure, or None if data is invalid.
    """
    if data is None:
        return None

    times = [datetime.fromisoformat(t) for t in data['times_iso']]
    n_times = len(times)
    blocks = data['blocks']

    fig = make_subplots(rows=min(len(blocks), 4), cols=len(blocks) // 4 + 1,
                        subplot_titles=[f"Elevations for {sb}" for sb in blocks] if len(blocks) > 1 else '')

    # (row, col, ant_names) per subplot so each one gets its own y-axis labels
    subplot_ant_names: list[tuple[int, int, list[str]]] = []
    for src_i, (src_block, bdata) in enumerate(blocks.items()):
        ant_names = bdata['ant_names']
        n_ants = len(ant_names)
        subplot_ant_names.append((src_i % 4 + 1, src_i // 4 + 1, ant_names))
        z_matrix = np.full((n_ants, n_times), np.nan)
        for anti, ant in enumerate(ant_names):
            vis = np.array(bdata['observability'][ant], dtype=bool)
            elev = np.array(bdata['elevations'][ant])
            z_matrix[n_ants - 1 - anti, vis] = elev[vis]

        station_info = data['station_info']
        hover_text = [[f"<b>{station_info[ant_names[n_ants-1-ai]]['name']}</b><br>"
                       f"Elevation: {z_matrix[ai, ti]:.0f}º<br>"
                       f"Time: {times[ti].strftime('%H:%M')}"
                       if not np.isnan(z_matrix[ai, ti]) else ""
                       for ti in range(n_times)] for ai in range(n_ants)]

        fig.add_trace(go.Heatmap(x=times, y=list(range(1, n_ants + 1)), z=z_matrix,
                                  colorscale='Viridis', zmin=5, zmax=90, showscale=show_colorbar,
                                  hoverinfo='text', text=hover_text, xgap=0, ygap=1),
                      row=src_i % 4 + 1, col=src_i // 4 + 1)

    axis_cfg = _get_axis_config_from_data(data)
    # Only attach the GST top axis when there is a single subplot, to avoid colliding with
    # the xaxis2/3/... slots that ``make_subplots`` reserves for additional source-block panels.
    if len(blocks) != 1:
        axis_cfg = {**axis_cfg, 'xaxis2_config': None}
    bottom_mirror = 'allticks' if axis_cfg.get('xaxis2_config') is None else False
    fig.update_layout(
        showlegend=False, hovermode='closest', xaxis_title=axis_cfg['xaxis_title'],
        yaxis_title="Antennas", title='', paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(type="date", tickformat="%H:%M", showline=True, linecolor='black', linewidth=1,
                   mirror=bottom_mirror, ticks='inside', tickmode='auto', range=axis_cfg['xaxis_range'],
                   showgrid=False, zeroline=False,
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        margin=dict(l=2, r=2, t=45, b=0))
    # Style all y axes, then label each subplot with its own block's antenna names.
    fig.update_yaxes(showline=True, linecolor='black', linewidth=1, mirror='allticks', ticks='inside',
                     showgrid=False, zeroline=False,
                     minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False))
    for row, col, names in subplot_ant_names:
        fig.update_yaxes(tickmode='array', tickvals=list(range(1, len(names) + 1)),
                         ticktext=list(reversed(names)), row=row, col=col)
    _apply_gst_top_axis(fig, axis_cfg)
    return fig


def elevation_curves_from_data(data: dict) -> Optional[go.Figure]:
    """Render the elevation curves plot from serialized data.

    Parameters
    ----------
    data : dict
        Serialized elevation data from serialize_elevation_data.

    Returns
    -------
    go.Figure or None
        Plotly figure, or None if data is invalid.
    """
    if data is None:
        return None

    times = [datetime.fromisoformat(t) for t in data['times_iso']]
    n_times = len(times)
    blocks = data['blocks']

    fig = make_subplots()
    fig.add_shape(type="rect", xref="paper", yref="y", x0=0, y0=0, x1=1, y1=20,
                  fillcolor="gray", opacity=0.2, layer="below", line_width=0)
    fig.add_shape(type="rect", xref="paper", yref="y", x0=0, y0=0, x1=1, y1=10,
                  fillcolor="gray", opacity=0.2, layer="below", line_width=0)

    station_info = data['station_info']
    for src_block, bdata in blocks.items():
        for ant in bdata['ant_names']:
            vis = np.array(bdata['observability'][ant], dtype=bool)
            y = np.full(n_times, None, dtype=float)
            y[vis] = np.array(bdata['elevations'][ant])[vis]
            sinfo = station_info[ant]
            fig.add_trace(go.Scatter(
                x=times, y=y, mode='lines', connectgaps=False, hoverinfo='none',
                hovertemplate=f"<b>{sinfo['name']} ({sinfo['codename']})</b><br><b>Time</b>: %{{x}}",
                name=sinfo['name']))

    axis_cfg = _get_axis_config_from_data(data)
    bottom_mirror = 'allticks' if axis_cfg.get('xaxis2_config') is None else False
    fig.update_layout(
        showlegend=True, hovermode='closest', xaxis_title=axis_cfg['xaxis_title'],
        yaxis_title="Elevation (degrees)", title='', paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(type="date", tickformat="%H:%M", showline=True, linecolor='black', linewidth=1,
                   mirror=bottom_mirror, ticks='inside', tickmode='auto', range=axis_cfg['xaxis_range'],
                   showgrid=False, zeroline=False,
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        yaxis=dict(showline=True, linecolor='black', linewidth=1, mirror='allticks', ticks='inside',
                   tickmode='auto', range=[0, 90], showgrid=False, zeroline=False,
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        margin=dict(l=2, r=2, t=45, b=0),
        legend=dict(x=0.01, y=0.99, xanchor='left', yanchor='top', bgcolor='rgba(255, 255, 255, 0.7)',
                    bordercolor='rgba(0, 0, 0, 0.3)', borderwidth=1))
    _apply_gst_top_axis(fig, axis_cfg)
    return fig


def serialize_worldmap_data(o) -> Optional[dict]:
    """Serialize station location/observability data for deferred worldmap rendering.

    Returns a JSON-serializable dict with all data needed by worldmap_from_data.

    Parameters
    ----------
    o : VLBIObs
        VLBI observation object.

    Returns
    -------
    dict or None
        Serialized worldmap data, or None if observation is invalid.
    """
    if o is None:
        return None

    try:
        ant_observes = o.can_be_observed()[list(o.can_be_observed().keys())[0]]
    except (ValueError, IndexError):
        ant_observes = {ant.codename: True for ant in o.stations}

    stations = []
    for ant in o.stations:
        stations.append({
            'lat': float(ant.location.lat.value),
            'lon': float(ant.location.lon.value),
            'name': ant.name,
            'country': ant.country,
            'diameter': ant.diameter,
            'observes': bool(ant_observes[ant.codename]),
        })
    return {'stations': stations}


def worldmap_from_data(data: dict) -> Optional[go.Figure]:
    """Render the worldmap from serialized data.

    Parameters
    ----------
    data : dict
        Serialized worldmap data from serialize_worldmap_data.

    Returns
    -------
    go.Figure or None
        Plotly figure, or None if data is invalid.
    """
    if data is None:
        return None

    stations = data['stations']
    lats = [s['lat'] for s in stations]
    lons = [s['lon'] for s in stations]
    colors = ['#a01d26' if s['observes'] else '#EAB308' for s in stations]
    hovertemplates = [f"{s['name']}<br>({s['country']})<br> {s['diameter']}<extra></extra>" for s in stations]
    avg_lon = np.mean(lons)

    fig = go.Figure(go.Scattergeo(lon=lons, lat=lats, mode='markers',
                                   marker=dict(size=10, color=colors), hovertemplate=hovertemplates))
    fig.update_geos(projection_type='orthographic', showland=True, landcolor='#9DB7C4',
                    projection_rotation=dict(lon=avg_lon, lat=0), bgcolor='rgba(0,0,0,0)')
    fig.update_layout(autosize=True, hovermode='closest', showlegend=False,
                      margin={'l': 0, 't': 0, 'b': 0, 'r': 0},
                      paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')
    return fig


def plot_worldmap_stations(o) -> Optional[go.Figure]:
    """Create a worldmap showing station locations and observability.

    Parameters
    ----------
    o : VLBIObs
        VLBI observation object.

    Returns
    -------
    go.Figure or None
        Plotly figure, or None if observation is invalid.
    """
    if o is None:
        return None

    data: dict[str, list] = defaultdict(list)
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
    avg_lon = np.mean(data['lon'])
    return _figure([dict(type='scattergeo', lon=data['lon'], lat=data['lat'], mode='markers',
                         marker=dict(size=10, color=['#a01d26' if q else '#EAB308' for q in data['observes']]),
                         # hover_name=data["text"], hover_data=None,
                         hovertemplate=data["hovertemplate"])],
                   dict(geo=dict(projection=dict(type='orthographic', rotation=dict(lon=avg_lon, lat=0)),
                                 showland=True, landcolor='#9DB7C4', bgcolor=_TRANSPARENT),
                        autosize=True, hovermode='closest', showlegend=False,
                        margin={'l': 0, 't': 0, 'b': 0, 'r': 0},
                        paper_bgcolor=_TRANSPARENT, plot_bgcolor=_TRANSPARENT))
