from collections import defaultdict
from typing import Optional
from datetime import datetime, timedelta
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from vlbiplanobs import sources


def _get_axis_config(o):
    """Returns axis configuration for elevation plots based on observation time mode.

    Inputs
    - o : VLBIObs
        VLBI observation object with times and gstimes attributes.

    Returns
    - dict with keys:
        - xaxis_range: [start, end] datetime range for x-axis
        - xaxis_title: title for x-axis
        - xaxis2_config: dict for secondary x-axis (or None if not needed)
    """
    time_range = [o.times.datetime[0], o.times.datetime[-1]]

    if o.fixed_time:
        # For xaxis (UTC): datetime values are used directly in time_range
        
        # For xaxis2 (GST): convert gstimes to datetime format for the range
        gst_start = datetime.combine(o.times.datetime[0].date(), datetime.min.time()) + timedelta(hours=o.gstimes[0].hour)
        gst_end = datetime.combine(o.times.datetime[0].date(), datetime.min.time()) + timedelta(hours=o.gstimes[-1].hour)
        gst_range = [gst_start, gst_end]
        
        # xaxis2_config = dict(overlaying='x', side='top', type='date', tickmode='auto',
        #                      showline=True, linecolor='black', linewidth=2, ticks='inside',
        #                      range=gst_range, title='Time (GST)', showgrid=False, zeroline=False,
        #                      tickformat='%H:%M')
        return {
            'xaxis_range': time_range,
            'xaxis_title': 'Time (UTC)',
            # 'xaxis2_range': gst_range,
            # 'xaxis2_config': xaxis2_config,
            # 'show_xaxis2': True
        }
    else:
        # No epoch defined: x-axis shows GST directly, x-axis2 shows GST too
        return {'xaxis_range': time_range, 'xaxis_title': 'Time (GST)'}


def elevation_plot_curves(o) -> Optional[go.Figure]:
    """Creates the plot showing when the different antennas can observe a given source,
    with the old style: with the elevation in the y-axis.
    """
    if o is None or not o.scans:
        return None

    srcup = o.is_observable()
    elevs = o.elevations()
    # fig = make_subplots(rows=min([len(srcup), 4]), cols=len(srcup) // 4 + 1,
    #                     subplot_titles=[f"Elevations for {src_block}" for src_block in srcup]
    #                     if len(srcup) > 1 else '')

    fig = make_subplots()

    fig.add_shape(type="rect",
                  xref="paper", yref="y",
                  x0=0, y0=0, x1=1, y1=20,
                  fillcolor="gray", opacity=0.2, layer="below", line_width=0)
    fig.add_shape(type="rect",
                  xref="paper", yref="y",
                  x0=0, y0=0, x1=1, y1=10,
                  fillcolor="gray", opacity=0.2, layer="below", line_width=0)
    for src_i, src_block in enumerate(srcup):
        for anti, ant in enumerate(srcup[src_block]):
            targets = o.scans[src_block].sources()  # sources.SourceType.TARGET)
            y = np.full(len(o.times.datetime), None, dtype=float)
            y[srcup[src_block][ant]] = elevs[targets[0].name][ant][srcup[src_block][ant]].value
            fig.add_trace(go.Scatter(x=o.times.datetime, y=y, mode='lines', connectgaps=False, hoverinfo='none',
                                     hovertemplate=f"<b>{o.stations[ant].name} ({o.stations[ant].codename})</b><br>"
                                                   "<b>Time</b>: %{x}",  # .strftime('%H:%M')}",
                                     name=o.stations[ant].name))

    # Get axis configuration based on observation time mode
    axis_cfg = _get_axis_config(o)

    # Add invisible trace for secondary x-axis (GST) only when needed
    # if axis_cfg['show_xaxis2']:
    #     fig.add_trace(go.Scatter(x=o.times.datetime, y=[None]*len(o.times.datetime),
    #                              mode='markers', marker=dict(opacity=0), showlegend=False,
    #                              xaxis='x2', hoverinfo='skip'))

    layout_kwargs = dict(showlegend=True, hovermode='closest', xaxis_title=axis_cfg['xaxis_title'],
                         yaxis_title="Elevation (degrees)", title='', paper_bgcolor='rgba(0,0,0,0)',
                         plot_bgcolor='rgba(0,0,0,0)',
                         xaxis=dict(type="date", tickformat="%H:%M", showline=True, linecolor='black', linewidth=1,
                                    mirror='allticks', ticks='inside', tickmode='auto', range=axis_cfg['xaxis_range'],
                                    showgrid=False, zeroline=False,
                                    minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
                         yaxis=dict(showline=True, linecolor='black', linewidth=1, mirror='allticks',
                                    ticks='inside', tickmode='auto', range=[0, 90],
                                    showgrid=False, zeroline=False,
                                    minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
                         margin=dict(l=2, r=2, t=45, b=0),
                         legend=dict(x=0.01, y=0.99, xanchor='left', yanchor='top',
                                     bgcolor='rgba(255, 255, 255, 0.7)',
                                     bordercolor='rgba(0, 0, 0, 0.3)', borderwidth=1))

    fig.update_layout(**layout_kwargs)
    return fig


def elevation_plot(o, show_colorbar: bool = False) -> Optional[go.Figure]:
    """Creates the plot showing when the different antennas can observe a given source.
    Optimized version using Heatmap instead of individual traces.
    """
    if o is None or not o.scans:
        return None

    srcup = o.is_observable()
    elevs = o.elevations()

    fig = make_subplots(rows=min([len(srcup), 4]), cols=len(srcup) // 4 + 1,
                        subplot_titles=[f"Elevations for {src_block}" for src_block in srcup]
                        if len(srcup) > 1 else '')

    for src_i, src_block in enumerate(srcup):
        ant_names = list(srcup[src_block].keys())
        n_ants = len(ant_names)
        n_times = len(o.times.datetime)

        # Build elevation matrix (antennas x times)
        z_matrix = np.full((n_ants, n_times), np.nan)
        for anti, ant in enumerate(ant_names):
            targets = o.scans[src_block].sources(sources.SourceType.TARGET)
            if len(targets) > 0:
                elev_values = elevs[targets[0].name][ant].value
            else:
                elev_values = elevs[o.scans[src_block].sources()[0].name][ant].value

            visibility = np.array(srcup[src_block][ant], dtype=bool)
            z_matrix[n_ants - 1 - anti, visibility] = elev_values[visibility]

        # Create hover text matrix
        hover_text = [[f"<b>{o.stations[ant_names[n_ants - 1 - ai]].name}</b><br>"
                       f"Elevation: {z_matrix[ai, ti]:.0f}º<br>"
                       f"Time: {o.times.datetime[ti].strftime('%H:%M')}"
                       if not np.isnan(z_matrix[ai, ti]) else ""
                       for ti in range(n_times)] for ai in range(n_ants)]

        fig.add_trace(
            go.Heatmap(
                x=o.times.datetime,
                y=list(range(1, n_ants + 1)),
                z=z_matrix,
                colorscale='Viridis',
                zmin=5, zmax=90,
                showscale=show_colorbar,
                hoverinfo='text',
                text=hover_text,
                xgap=0, ygap=1
            ),
            row=src_i % 4 + 1,
            col=src_i // 4 + 1)

    # Get axis configuration based on observation time mode
    axis_cfg = _get_axis_config(o)

    layout_kwargs = dict(
        showlegend=False,
        hovermode='closest',
        xaxis_title=axis_cfg['xaxis_title'],
        yaxis_title="Antennas",
        title='',
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(type="date", tickformat="%H:%M", showline=True, linecolor='black', linewidth=1,
                   mirror='allticks', ticks='inside', tickmode='auto', range=axis_cfg['xaxis_range'],
                   showgrid=False, zeroline=False,
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        yaxis=dict(showline=True, linecolor='black', linewidth=1, mirror='allticks', ticks='inside',
                   tickmode='array',
                   tickvals=list(range(1, n_ants + 1)),
                   ticktext=list(reversed(ant_names)),
                   showgrid=False, zeroline=False,
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        margin=dict(l=2, r=2, t=45, b=0))

    fig.update_layout(**layout_kwargs)
    if show_colorbar:
        fig.update_layout(coloraxis=dict(colorscale='Viridis'),
                          coloraxis_colorbar=dict(title='Elevation (degrees)'))

    return fig


def serialize_uv_data(o) -> Optional[dict]:
    """Serialize UV data for storage in dcc.Store.
    Returns dict of baseline -> {'x': [...], 'y': [...]} with both +/- uv points.
    """
    if o is None or not o.scans:
        return None
    
    bl_uv = o.get_uv_data()
    serialized = {}
    for key, arr in bl_uv[list(bl_uv.keys())[0]].items():
        arr_values = arr.value if hasattr(arr, 'value') else arr
        serialized[key] = {
            'x': arr_values[:, 0].tolist() + (-arr_values[:, 0]).tolist(),
            'y': arr_values[:, 1].tolist() + (-arr_values[:, 1]).tolist()
        }
    return serialized


def uvplot_from_data(uv_data: dict, filter_antennas: Optional[list[str]] = None) -> Optional[go.Figure]:
    """Creates UV coverage plot from serialized UV data."""
    if uv_data is None:
        return None
    
    highlight_colors = [
        "#FF0000", "#0000FF", "#008000", "#FFA500", "#800080",
        "#00FFFF", "#FF00FF", "#FFD700", "#00FF00", "#A52A2A",
        "#FFC0CB", "#808000", "#008080", "#000080", "#FF7F50",
        "#4B0082", "#FF8C00", "#40E0D0", "#6A5ACD", "#006400",
    ]

    # Group points by color for batching
    color_groups: dict[str, dict] = {'black': {'x': [], 'y': [], 'text': []}}
    if filter_antennas:
        for i, ant in enumerate(filter_antennas):
            color_groups[highlight_colors[i % len(highlight_colors)]] = {'x': [], 'y': [], 'text': []}

    def get_color(baseline: str) -> str:
        if filter_antennas:
            for i, ant in enumerate(filter_antennas):
                if ant in baseline:
                    return highlight_colors[i % len(highlight_colors)]
        return 'black'

    for baseline, points in uv_data.items():
        color = get_color(baseline)
        color_groups[color]['x'].extend(points['x'])
        color_groups[color]['y'].extend(points['y'])
        color_groups[color]['text'].extend([baseline] * len(points['x']))

    # Create one trace per color
    data = []
    for color, points in color_groups.items():
        if points['x']:
            data.append(go.Scattergl(
                x=points['x'],
                y=points['y'],
                mode='markers',
                marker=dict(color=color, size=2 if color == 'black' else 4),
                hovertext=points['text'],
                hoverinfo='text',
                showlegend=False
            ))

    fig = go.Figure(data=data)
    fig.update_layout(
        showlegend=False,
        hovermode='closest',
        xaxis_title='u (λ)',
        yaxis_title='v (λ)',
        title='',
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(showline=True, linecolor='black', linewidth=1, mirror='allticks', ticks='inside',
                   tickmode='auto', showgrid=False, zeroline=False,
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        yaxis=dict(showline=True, linecolor='black', linewidth=1, mirror='allticks', ticks='inside',
                   scaleanchor="x", scaleratio=1, tickmode='auto', showgrid=False, zeroline=False,
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        margin=dict(l=0, r=0, t=0, b=0)
    )
    fig.update_xaxes(constrain='domain')
    return fig


def uvplot(o, filter_antennas: Optional[list[str]] = None) -> Optional[go.Figure]:
    """Creates UV coverage plot. Optimized to batch all points into fewer traces."""
    if o is None or not o.scans:
        return None

    bl_uv = o.get_uv_data()
    highlight_colors = [
        "#FF0000", "#0000FF", "#008000", "#FFA500", "#800080",
        "#00FFFF", "#FF00FF", "#FFD700", "#00FF00", "#A52A2A",
        "#FFC0CB", "#808000", "#008080", "#000080", "#FF7F50",
        "#4B0082", "#FF8C00", "#40E0D0", "#6A5ACD", "#006400",
    ]

    # Group points by color for batching
    color_groups: dict[str, dict] = {'black': {'x': [], 'y': [], 'text': []}}
    if filter_antennas:
        for i, ant in enumerate(filter_antennas):
            color_groups[highlight_colors[i % len(highlight_colors)]] = {'x': [], 'y': [], 'text': []}

    def get_color(baseline: str) -> str:
        if filter_antennas:
            for i, ant in enumerate(filter_antennas):
                if ant in baseline:
                    return highlight_colors[i % len(highlight_colors)]
        return 'black'

    for key, arr in bl_uv[list(bl_uv.keys())[0]].items():
        color = get_color(key)
        # Add both +uv and -uv points (use .value for Quantity arrays)
        arr_values = arr.value if hasattr(arr, 'value') else arr
        color_groups[color]['x'].extend(arr_values[:, 0].tolist())
        color_groups[color]['x'].extend((-arr_values[:, 0]).tolist())
        color_groups[color]['y'].extend(arr_values[:, 1].tolist())
        color_groups[color]['y'].extend((-arr_values[:, 1]).tolist())
        color_groups[color]['text'].extend([key] * (2 * len(arr_values)))

    # Create one trace per color (much fewer traces)
    data = []
    for color, points in color_groups.items():
        if points['x']:
            data.append(go.Scattergl(  # Use Scattergl for WebGL rendering - faster
                x=points['x'],
                y=points['y'],
                mode='markers',
                marker=dict(color=color, size=2 if color == 'black' else 4),
                hovertext=points['text'],
                hoverinfo='text',
                showlegend=False
            ))

    fig = go.Figure(data=data)
    fig.update_layout(
        showlegend=False,
        hovermode='closest',
        xaxis_title='u (λ)',
        yaxis_title='v (λ)',
        title='',
        paper_bgcolor='rgba(0,0,0,0)',   # Transparent background
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(showline=True, linecolor='black', linewidth=1, mirror='allticks', ticks='inside',
                   tickmode='auto', showgrid=False, zeroline=False,
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        yaxis=dict(showline=True, linecolor='black', linewidth=1, mirror='allticks', ticks='inside',
                   scaleanchor="x", scaleratio=1, tickmode='auto', showgrid=False, zeroline=False,  # Square aspect ratio
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        margin=dict(l=0, r=0, t=0, b=0)
    )
    fig.update_xaxes(constrain='domain')
    return fig


def serialize_elevation_data(o) -> Optional[dict]:
    """Serialize elevation/observability data for deferred plot rendering.

    Returns a JSON-serializable dict with all data needed by elevation_plot_from_data
    and elevation_curves_from_data. Returns None if no scans.
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

        observability = {}
        elevation_vals = {}
        for ant in ant_names:
            observability[ant] = [bool(v) for v in srcup[src_block][ant]]
            elev_arr = elevs[target_name][ant].value
            elevation_vals[ant] = [float(v) for v in elev_arr]

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
    """Build axis config from serialized elevation data."""
    times = [datetime.fromisoformat(t) for t in data['times_iso']]
    time_range = [times[0], times[-1]]
    if data['fixed_time']:
        return {'xaxis_range': time_range, 'xaxis_title': 'Time (UTC)'}
    else:
        return {'xaxis_range': time_range, 'xaxis_title': 'Time (GST)'}


def elevation_plot_from_data(data: dict, show_colorbar: bool = False) -> Optional[go.Figure]:
    """Render the heatmap elevation plot from serialized data (no obs object needed)."""
    if data is None:
        return None

    times = [datetime.fromisoformat(t) for t in data['times_iso']]
    n_times = len(times)
    blocks = data['blocks']

    fig = make_subplots(rows=min(len(blocks), 4), cols=len(blocks) // 4 + 1,
                        subplot_titles=[f"Elevations for {sb}" for sb in blocks] if len(blocks) > 1 else '')

    for src_i, (src_block, bdata) in enumerate(blocks.items()):
        ant_names = bdata['ant_names']
        n_ants = len(ant_names)
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
    fig.update_layout(
        showlegend=False, hovermode='closest', xaxis_title=axis_cfg['xaxis_title'],
        yaxis_title="Antennas", title='', paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(type="date", tickformat="%H:%M", showline=True, linecolor='black', linewidth=1,
                   mirror='allticks', ticks='inside', tickmode='auto', range=axis_cfg['xaxis_range'],
                   showgrid=False, zeroline=False,
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        yaxis=dict(showline=True, linecolor='black', linewidth=1, mirror='allticks', ticks='inside',
                   tickmode='array', tickvals=list(range(1, n_ants + 1)), ticktext=list(reversed(ant_names)),
                   showgrid=False, zeroline=False,
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        margin=dict(l=2, r=2, t=45, b=0))
    return fig


def elevation_curves_from_data(data: dict) -> Optional[go.Figure]:
    """Render the elevation curves plot from serialized data (no obs object needed)."""
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
    fig.update_layout(
        showlegend=True, hovermode='closest', xaxis_title=axis_cfg['xaxis_title'],
        yaxis_title="Elevation (degrees)", title='', paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(type="date", tickformat="%H:%M", showline=True, linecolor='black', linewidth=1,
                   mirror='allticks', ticks='inside', tickmode='auto', range=axis_cfg['xaxis_range'],
                   showgrid=False, zeroline=False,
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        yaxis=dict(showline=True, linecolor='black', linewidth=1, mirror='allticks', ticks='inside',
                   tickmode='auto', range=[0, 90], showgrid=False, zeroline=False,
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        margin=dict(l=2, r=2, t=45, b=0),
        legend=dict(x=0.01, y=0.99, xanchor='left', yanchor='top', bgcolor='rgba(255, 255, 255, 0.7)',
                    bordercolor='rgba(0, 0, 0, 0.3)', borderwidth=1))
    return fig


def serialize_worldmap_data(o) -> Optional[dict]:
    """Serialize station location/observability data for deferred worldmap rendering.

    Returns a JSON-serializable dict with all data needed by worldmap_from_data.
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
    """Render the worldmap from serialized data (no obs object needed)."""
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
    fig = go.Figure(go.Scattergeo(
        lon=data['lon'],
        lat=data['lat'],
        mode='markers',
        marker=dict(size=10, color=['#a01d26' if q else '#EAB308' for q in data['observes']]),
        # hover_name=data["text"], hover_data=None,
        hovertemplate=data["hovertemplate"])
    )

    fig.update_geos(projection_type='orthographic', showland=True, landcolor='#9DB7C4',
                    projection_rotation=dict(lon=avg_lon, lat=0), bgcolor='rgba(0,0,0,0)')

    fig.update_layout(autosize=True, hovermode='closest', showlegend=False,
                      margin={'l': 0, 't': 0, 'b': 0, 'r': 0},
                      paper_bgcolor='rgba(0,0,0,0)',
                      plot_bgcolor='rgba(0,0,0,0)')
    return fig
