import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import plotly.graph_objects as go

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
    assert o.times is not None
    # Normalize the color array to [0, 1] and map it to the viridis colormap
    norm = Normalize(vmin=0, vmax=90)  # Normalize values between 0 and 90
    viridis = cm.get_cmap('viridis')

    # Create figure
    fig: dict[str, go.Figure] = {}

    srcup = o.is_observable()
    elevs = o.elevations()
    # Break the line into segments for color changes
    for src_block in srcup:
        fig[src_block] = go.Figure()
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

            y_value = np.ones(2) * anti
            for i in range(len(srcup[src_block][ant][srcup[src_block][ant]]) - 1):
                fig[src_block].add_trace(
                    go.Scatter(
                        x=o.times.datetime[srcup[src_block][ant]][i:i+2],
                        y=y_value,
                        mode="lines",
                        line=dict(color=color_str[i], width=10),
                        showlegend=False,
                        hovertemplate=f"Elevation: {colors[i]:.0f}º<extra></extra>"
                    )
                )

        # Update layout
        fig[src_block].update_layout(
            xaxis_title="Time",
            yaxis_title="Antennas",
            xaxis=dict(type="date"),
            title=f"Elevations for {src_block}"
        )
        fig[src_block].update_yaxes(
            tickvals=list(range(len(srcup[src_block].keys()))),  # Positions on the y-axis
            ticktext=list(srcup[src_block].keys()),   # Corresponding labels (names)
            title="Antennas"
        )

    return fig
