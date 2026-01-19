from collections import defaultdict
from typing import Optional
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from vlbiplanobs import sources


def _gst_tickvals_ticktext(o):
    """Returns tickvals and ticktext for a GST secondary x-axis.
    
    The tickvals are datetime positions on the x-axis, and ticktext shows
    the corresponding GST time at each position.
    
    When fixed_time is False, the x-axis datetimes already represent GST hours
    (synthetic times), so both axes should show the same values.
    When fixed_time is True, the x-axis shows real UTC times, and we need to
    show the corresponding GST times on the secondary axis.
    """
    n_times = len(o.times.datetime)
    
    # Select evenly spaced indices across the full time range
    n_ticks = min(8, n_times)
    if n_ticks <= 1:
        indices = [0]
    else:
        indices = [int(i * (n_times - 1) / (n_ticks - 1)) for i in range(n_ticks)]
    
    # tickvals are the x-axis positions (datetime)
    tickvals = [o.times.datetime[i] for i in indices]
    
    if o.fixed_time:
        # Real UTC times - show corresponding GST
        gst = o.gstimes.hour
        ticktext = [f"{int(gst[i]):02d}:{int((gst[i] % 1) * 60):02d}" for i in indices]
    else:
        # Synthetic times representing GST - extract hour:minute from datetime
        ticktext = [o.times.datetime[i].strftime('%H:%M') for i in indices]
    
    return tickvals, ticktext


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
            fig.add_trace(
                go.Scatter(x=o.times.datetime,
                           y=y,
                           mode='lines',
                           connectgaps=False,
                           hoverinfo='none',
                           hovertemplate=f"<b>{o.stations[ant].name} ({o.stations[ant].codename})</b><br>"
                                         "<b>Time</b>: %{x}",  # .strftime('%H:%M')}",
                           name=o.stations[ant].name))

    # Add invisible trace for secondary x-axis (GST) to make it visible
    tickvals, ticktext = _gst_tickvals_ticktext(o)
    fig.add_trace(go.Scatter(x=o.times.datetime, y=[None]*len(o.times.datetime),
                             mode='markers', marker=dict(opacity=0), showlegend=False,
                             xaxis='x2', hoverinfo='skip'))

    layout_kwargs = dict(
        showlegend=True,
        hovermode='closest',
        xaxis_title="Time (GST)" if not o.fixed_time else "Time (UTC)",
        yaxis_title="Elevation (degrees)",
        title='',
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(type="date", tickformat="%H:%M", showline=True, linecolor='black', linewidth=1,
                   mirror=False, ticks='inside', tickmode='auto',
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        xaxis2=dict(overlaying='x', side='top', type='date', tickmode='array',
                    tickvals=tickvals, ticktext=ticktext, showline=True,
                    linecolor='black', linewidth=1, ticks='inside',
                    title='Time (GST)'),
        yaxis=dict(showline=True, linecolor='black', linewidth=1, mirror='allticks', ticks='inside',
                   tickmode='auto', range=[0, 90],
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        margin=dict(l=2, r=2, t=45, b=0),
        legend=dict(x=0.01, y=0.99, xanchor='left', yanchor='top', bgcolor='rgba(255, 255, 255, 0.7)',
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

    # Add invisible trace for secondary x-axis (GST) to make it visible
    tickvals, ticktext = _gst_tickvals_ticktext(o)
    fig.add_trace(go.Scatter(x=o.times.datetime, y=[None]*len(o.times.datetime),
                             mode='markers', marker=dict(opacity=0), showlegend=False,
                             xaxis='x2', hoverinfo='skip'))

    layout_kwargs = dict(
        showlegend=False,
        hovermode='closest',
        xaxis_title="Time (UTC)" if o.fixed_time else "Time (GST)",
        yaxis_title="Antennas",
        title='',
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(type="date", tickformat="%H:%M", showline=True, linecolor='black', linewidth=1,
                   mirror=False, ticks='inside', tickmode='auto',
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        xaxis2=dict(overlaying='x', side='top', type='date', tickmode='array',
                    tickvals=tickvals, ticktext=ticktext, showline=True,
                    linecolor='black', linewidth=1, mirror=False, ticks='inside',
                    minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False),
                    title='Time (GST)'),
        yaxis=dict(showline=True, linecolor='black', linewidth=1, mirror='allticks', ticks='inside',
                   tickmode='array',
                   tickvals=list(range(1, n_ants + 1)),
                   ticktext=list(reversed(ant_names)),
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
        xaxis=dict(
            showline=True,
            linecolor='black',
            linewidth=1,
            mirror='allticks',
            ticks='inside',
            tickmode='auto',
            minor=dict(
                ticks='inside',
                ticklen=4,
                tickcolor='black',
                showgrid=False
            )
        ),
        yaxis=dict(
            showline=True,
            linecolor='black',
            linewidth=1,
            mirror='allticks',
            ticks='inside',
            scaleanchor="x",
            scaleratio=1,
            tickmode='auto',
            minor=dict(
                ticks='inside',
                ticklen=4,
                tickcolor='black',
                showgrid=False
            )
        ),
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
        xaxis=dict(
            showline=True,
            linecolor='black',
            linewidth=1,
            mirror='allticks',
            ticks='inside',
            tickmode='auto',
            minor=dict(
                ticks='inside',
                ticklen=4,
                tickcolor='black',
                showgrid=False
            )
        ),
        yaxis=dict(
            showline=True,
            linecolor='black',
            linewidth=1,
            mirror='allticks',
            ticks='inside',
            scaleanchor="x",      # Square aspect ratio
            scaleratio=1,
            tickmode='auto',
            minor=dict(
                ticks='inside',
                ticklen=4,
                tickcolor='black',
                showgrid=False
            )
        ),
        margin=dict(l=0, r=0, t=0, b=0)
    )
    fig.update_xaxes(constrain='domain')
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

    fig.update_geos(
        # projection_type='natural earth',
        projection_type='orthographic',
        # lonaxis_range=[avg_lon-90, avg_lon+90],
        showland=True,
        landcolor='#9DB7C4',
        projection_rotation=dict(lon=avg_lon, lat=0),
        bgcolor='rgba(0,0,0,0)'
    )

    fig.update_layout(autosize=True, hovermode='closest', showlegend=False,
                      margin={'l': 0, 't': 0, 'b': 0, 'r': 0},
                      paper_bgcolor='rgba(0,0,0,0)',
                      plot_bgcolor='rgba(0,0,0,0)')
    return fig
