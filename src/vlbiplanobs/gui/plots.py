import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from astropy.time import Time
from astropy import units as u
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from vlbiplanobs import observation as obs
from vlbiplanobs import sources

# Example data
data = {
    "time": pd.date_range(start="2025-03-01", periods=10, freq="D"),
    "name": ["Alice"] * 10,
    "value": np.random.rand(10),
    "color": np.linspace(0, 90, 10)  # Example color array with values between 0 and 90
}


def elevation_plot(o: obs.Observation):
    """Creates the plot showing when the different antennas can observe a given
    source,
    """
    assert o.scans is not None
    # Normalize the color array to [0, 1] and map it to the viridis colormap
    norm = Normalize(vmin=0, vmax=90)  # Normalize values between 0 and 90
    viridis = cm.get_cmap('viridis')

    # Create figure
    # fig: dict[str, go.Figure] = {}

    if o.times is None:
        # Let's calculate the visibility GST ranges
        localtimes = Time('2025-09-21', scale='utc') + np.arange(0.0, 1.005, 0.01)*u.day
        o.times = localtimes
        doing_gst = True
    else:
        doing_gst = False

    srcup = o.is_observable()
    elevs = o.elevations()
    if doing_gst:
        o.times = None

    fig = make_subplots(rows=min([len(srcup), 4]), cols=len(srcup) // 4 + 1,
                        subplot_titles=[f"Elevations for {src_block}" for src_block in srcup])
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
                # fig[src_block].add_trace(
                fig.add_trace(
                    go.Scatter(
                        x=o.times.datetime[srcup[src_block][ant]][i:i+2] if o.times is not None else
                                                  localtimes.datetime[srcup[src_block][ant]][i:i+2],
                        y=y_value,
                        mode="lines",
                        line=dict(color=color_str[i], width=10),
                        showlegend=False,
                        hovertemplate=f"<b>Elevation</b>: {colors[i]:.0f}º<extra></extra><br>"
                                      f"<b>Time</b>: {o.times[i].strftime('%H:%M') if o.times is not None
                                                      else localtimes[i].strftime('%H:%M')}",
                    ),
                    row=src_i % 4 + 1,
                    col=src_i // 4 + 1
                )

        # Update layout
        fig.update_layout(
            xaxis_title="Time (GST)" if doing_gst else "Time (UTC)",
            yaxis_title="Antennas",
            xaxis=dict(type="date", tickformat="%H:%M"),
            # title=f"Elevations for {src_block}"
        )
        fig.update_yaxes(
            tickvals=len(o.stations) - np.array(range(len(srcup[src_block].keys()))),  # Positions on the y-axis
            # tickvals=list(range(len(srcup[src_block].keys()))),  # Positions on the y-axis
            ticktext=list(srcup[src_block].keys()),   # Corresponding labels (names)
            title="Antennas"
        )

    return fig
