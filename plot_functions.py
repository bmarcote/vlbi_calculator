import numpy as np
import plotly.graph_objs as go




def time_elevation(times, elevations, ant_names, xlim=None):
    """
    """
    data = []
    for time,elev,ant in zip(times, elevations, ant_names):
        data.append(go.Scatter(x=time.datetime, y=elev, mode='markers', name=ant))

    layout = go.Layout(title='Source elevation vs time', showlegend=True,
                       xaxis={'title': 'Time (UTC)', 'range': xlim},
                       yaxis={'title': 'Elevation (deg)'}, hovermode='closest')
    return {'data': data, 'layout': layout}


def time_visibility(times, elevations, ant_names, xlim=None):
    data = []
    for i,time,elev,ant in zip(np.arange(len(ant_names)), times, elevations, ant_names):
        data.append(go.Scatter(x=time.datetime, y=np.zeros_like(elev).deg+i, mode='markers',
                               name=ant))

    layout = go.Layout(title='Visibility of source vs time', showlegend=False,
                       xaxis={'title': 'Time (UTC)', 'range': xlim},
                       yaxis={'title': 'Antennas', 'tickmode': 'array',
                              'ticktext': ant_names, 'tickvals': np.arange(len(ant_names))},
                       hovermode=False)
    return {'data': data, 'layout': layout}


