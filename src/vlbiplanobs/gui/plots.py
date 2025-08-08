from collections import defaultdict
from typing import Optional
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from vlbiplanobs import sources


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

    fig.update_layout(
        showlegend=True,
        # showgrid=False,
        hovermode='closest',
        xaxis_title="Time (GST)" if not o.fixed_time else "Time (UTC)",
        yaxis_title="Elevation (degrees)",
        title='',
        paper_bgcolor='rgba(0,0,0,0)',   # Transparent background
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(type="date", tickformat="%H:%M", showline=True, linecolor='black', linewidth=1,
                   mirror='allticks', ticks='inside', tickmode='auto',
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        yaxis=dict(showline=True, linecolor='black', linewidth=1, mirror='allticks', ticks='inside',
                   tickmode='auto', range=[0, 90],
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        margin=dict(l=2, r=2, t=1, b=0),
        legend=dict(x=0.01, y=0.99, xanchor='left', yanchor='top', bgcolor='rgba(255, 255, 255, 0.7)',
                    bordercolor='rgba(0, 0, 0, 0.3)', borderwidth=1))

    return fig


def elevation_plot(o, show_colorbar: bool = False) -> Optional[go.Figure]:
    """Creates the plot showing when the different antennas can observe a given
    source,
    """
    if o is None or not o.scans:
        return None

    # Normalize the color array to [0, 1] and map it to the viridis colormap
    norm = Normalize(vmin=5, vmax=90)  # Normalize values between 0 and 90
    viridis = cm.get_cmap('viridis')
    srcup = o.is_observable()
    elevs = o.elevations()

    fig = make_subplots(rows=min([len(srcup), 4]), cols=len(srcup) // 4 + 1,
                        subplot_titles=[f"Elevations for {src_block}" for src_block in srcup]
                        if len(srcup) > 1 else '')

    for src_i, src_block in enumerate(srcup):
        for anti, ant in enumerate(srcup[src_block]):
            targets = o.scans[src_block].sources(sources.SourceType.TARGET)
            if len(targets) > 0:
                colors = elevs[targets[0].name][ant].value
            else:
                colors = elevs[o.scans[src_block].sources()[0].name][ant].value

            colors_cm = viridis(norm(colors))
            color_str = np.array([f"rgba({r}, {g}, {b}, {a})" for r, g, b, a in colors_cm])
            visibility = np.array(srcup[src_block][ant], dtype=bool)
            color_str[~visibility] = 'rgba(0, 0, 0, 0)'
            y_value = np.zeros(2) + (len(o.stations) - anti)
            visible_indices = np.where(visibility[:-1] & visibility[1:])[0]
            for i in visible_indices:
                fig.add_trace(
                    go.Scatter(
                        x=o.times.datetime[i:i + 2],
                        y=y_value,
                        mode="lines",
                        line=dict(color=color_str[i], width=10),
                        showlegend=False,
                        marker=dict(showscale=show_colorbar),
                        hovertemplate=f"<b>{o.stations[ant].name} ({o.stations[ant].codename})</b><br>"
                                    "<b>Elevation</b>: "
                                    f"{colors[i]:.0f}º<extra></extra><br><b>Time</b>: "
                                    f"{o.times.datetime[i].strftime('%H:%M')}"),
                    row=src_i % 4 + 1,
                    col=src_i // 4 + 1)

    fig.update_layout(
        showlegend=False,
        hovermode='closest',
        xaxis_title="Time (GST)" if not o.fixed_time else "Time (UTC)",
        yaxis_title="Antennas",
        title='',
        paper_bgcolor='rgba(0,0,0,0)',   # Transparent background
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(type="date", tickformat="%H:%M", showline=True, linecolor='black', linewidth=1,
                   mirror='allticks', ticks='inside', tickmode='auto',
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        yaxis=dict(showline=True, linecolor='black', linewidth=1, mirror='allticks', ticks='inside',
                   tickmode='array',
                   tickvals=len(o.stations) - np.array(range(len(srcup[src_block].keys()))),
                   ticktext=list(srcup[src_block].keys()),
                   minor=dict(ticks='inside', ticklen=4, tickcolor='black', showgrid=False)),
        margin=dict(l=2, r=2, t=1, b=0))
    if show_colorbar:
        fig.update_layout(coloraxis=dict(colorscale='Viridis'),
                          coloraxis_colorbar=dict(title='Elevation (degrees)'))

    return fig


def uvplot(o, filter_antennas: Optional[list[str]] = None) -> Optional[go.Figure]:
    if o is None or not o.scans:
        return None

    bl_uv = o.get_uv_data()
    data = []
    highlight_colors = [
        "#FF0000",  # Red
        "#0000FF",  # Blue
        "#008000",  # Green
        "#FFA500",  # Orange
        "#800080",  # Purple
        "#00FFFF",  # Cyan
        "#FF00FF",  # Magenta
        "#FFD700",  # Gold
        "#00FF00",  # Lime
        "#A52A2A",  # Brown
        "#FFC0CB",  # Pink
        "#808000",  # Olive
        "#008080",  # Teal
        "#000080",  # Navy
        "#FF7F50",  # Coral
        "#4B0082",  # Indigo
        "#FF8C00",  # DarkOrange
        "#40E0D0",  # Turquoise
        "#6A5ACD",  # SlateBlue
        "#006400",  # DarkGreen
    ]
    if filter_antennas is not None and len(filter_antennas) > len(highlight_colors):
        highlight_colors = highlight_colors * (len(filter_antennas) // len(highlight_colors))

    def get_color(baseline: str, filter_antennas: list[str]):
        for i, ant in enumerate(filter_antennas):
            if ant in baseline:
                return highlight_colors[i]

        return 'black'

    for key, arr in bl_uv[list(bl_uv.keys())[0]].items():
        if filter_antennas is not None:
            color = get_color(key, filter_antennas)
        else:
            color = 'black'

        uv = np.empty((2*len(arr), 2))
        uv[:len(arr), :] = arr
        uv[len(arr):, :] = -arr
        hover_text = [key] * len(uv)
        data.append(go.Scatter(
            x=uv[:, 0],
            y=uv[:, 1],
            mode='markers',
            marker=dict(color=color, size=2 if color == 'black' else 4),
            name=key,
            hovertext=hover_text,
            hoverinfo='text'
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
        projection_rotation=dict(lon=avg_lon, lat=0)
    )

    fig.update_layout(autosize=True, hovermode='closest', showlegend=False,
                      margin={'l': 0, 't': 0, 'b': 0, 'r': 0})
    return fig
