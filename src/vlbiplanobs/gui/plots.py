from typing import Optional
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from astropy.time import Time
from astropy import units as u
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
from dash import html
from vlbiplanobs import observation as obs
from vlbiplanobs import sources


def elevation_plot(o, show_colorbar: bool = False) -> go.Figure:
    """Creates the plot showing when the different antennas can observe a given
    source,
    """
    if o is None:
        return html.Div()

    assert o.scans is not None
    # Normalize the color array to [0, 1] and map it to the viridis colormap
    norm = Normalize(vmin=0, vmax=90)  # Normalize values between 0 and 90
    viridis = cm.get_cmap('viridis')

    if o.times is None:
        localtimes = o._REF_TIMES
    else:
        localtimes = o.times

    srcup = o.is_observable_at(localtimes)
    elevs = o.elevations()  # auto does it in localtimes
    # if doing_gst:
    #     o.times = None

    fig = make_subplots(rows=min([len(srcup), 4]), cols=len(srcup) // 4 + 1,
                        subplot_titles=[f"Elevations for {src_block}" for src_block in srcup] \
                        if len(srcup) > 1 else '')

    # Break the line into segments for color changes
    for src_i, src_block in enumerate(srcup):
        # fig[src_block] = go.Figure()
        for anti, ant in enumerate(srcup[src_block]):
            targets = o.scans[src_block].sources(sources.SourceType.TARGET)
            if len(targets) > 0:
                colors = elevs[src_block][targets[0].name][ant][srcup[src_block][ant]].value
            else:
                colors = elevs[src_block][o.scans[src_block].sources()[0].name][ant][srcup[src_block][ant]].value

            colors_cm = viridis(norm(colors))

            # Convert RGBA colors to Plotly-compatible string format
            # color_str = colors #.apply(lambda c: f"rgba({int(c[0]*255)}, {int(c[1]*255)}, {int(c[2]*255)}, {c[3]})")
            color_str = [f"rgba({r}, {g}, {b}, {a})" for r, g, b, a in colors_cm]

            y_value = np.zeros(2) + (len(o.stations) - anti)
            for i in range(len(srcup[src_block][ant][srcup[src_block][ant]]) - 1):
                if colors[i] > 10:
                    # fig[src_block].add_trace(
                    fig.add_trace(
                        go.Scatter(
                            x=o.times.datetime[srcup[src_block][ant]][i:i+2] if o.times is not None else
                                                      localtimes.datetime[srcup[src_block][ant]][i:i+2],
                            y=y_value,
                            mode="lines",
                            line=dict(color=color_str[i], width=10),
                            showlegend=False,
                            marker=dict(showscale=show_colorbar),
                            hovertemplate=f"<b>{o.stations[ant].name}</b><br><b>Elevation</b>: {colors[i]:.0f}º<extra></extra><br>"
                                          f"<b>Time</b>: {o.times[i].strftime('%H:%M') \
                                          if o.times is not None \
                                                          else localtimes[i].strftime('%H:%M')}",
                        ),
                        row=src_i % 4 + 1,
                        col=src_i // 4 + 1
                    )

    fig.update_layout(
        showlegend=False,
        hovermode='closest',
        xaxis_title="Time (GST)" if o.times is None else "Time (UTC)",
        yaxis_title="Antennas",
        title='',
        paper_bgcolor='rgba(0,0,0,0)',   # Transparent background
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(
            type="date",
            tickformat="%H:%M",
            showline=True,
            linecolor='black',
            linewidth=1,
            mirror='allticks',  # True,
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
            tickmode='array',
            tickvals=len(o.stations) - np.array(range(len(srcup[src_block].keys()))),
            ticktext=list(srcup[src_block].keys()),
            minor=dict(
                ticks='inside',
                ticklen=4,
                tickcolor='black',
                showgrid=False
            )
        ),
        margin=dict(l=2, r=2, t=1, b=0)
    )
    if show_colorbar:
        fig.update_layout(coloraxis=dict(colorscale='Viridis'),
                          coloraxis_colorbar=dict(title='Elevation (degrees)'))

    return fig


def uvplot(o, filter_antennas: Optional[list[str]] = None) -> go.Figure:
    if o is None:
        return None

    assert o.scans is not None
    bl_uv = o.get_uv_data()
    data = []
    default_color = 'black'
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
            marker=dict(color=color, size=1 if color == 'black' else 2),
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


def plot_worldmap_stations(data):
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
        lonaxis_range=[avg_lon-90, avg_lon+90],
        showland=True,
        landcolor='#9DB7C4',
        center=dict(lon=avg_lon, lat=0)
    )

    fig.update_layout(autosize=True, hovermode='closest', showlegend=False,
                      margin={'l': 0, 't': 0, 'b': 0, 'r': 0})
    return fig
