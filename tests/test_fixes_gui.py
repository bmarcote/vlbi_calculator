"""Tests for the GUI server-side validation (gui/validation.py) and the uv-plot antenna highlight matching."""
import json
import math
import numpy as np
import plotly
import pytest
from vlbiplanobs import freqsetups as fs
from vlbiplanobs import observation
from vlbiplanobs.gui import validation as val
from vlbiplanobs.gui import plots


def _valid_params(**overrides) -> dict:
    """Return a valid store-obs-params-shaped dict, with optional overrides."""
    params = {'band': '18cm', 'stations': ['Ef', 'Mc', 'Nt'], 'targets': ['3C273'], 'duration': 10,
              'ontarget': 0.7, 'startdate': None, 'starttime': None, 'datarate': 2048, 'subbands': 8,
              'channels': 64, 'polarizations': 2, 'inttime': 2}
    params.update(overrides)
    return params


# ------------------------------------------------------------------------------------------
# Duration
# ------------------------------------------------------------------------------------------
def test_duration_bounds():
    assert val.validate_duration(None) is None
    assert val.validate_duration(0) is None
    assert val.validate_duration(-3) is None
    assert val.validate_duration(12) == 12.0
    assert val.validate_duration(val.MAX_DURATION_H) == val.MAX_DURATION_H
    with pytest.raises(val.InvalidObsParams):
        val.validate_duration(val.MAX_DURATION_H + 0.1)
    with pytest.raises(val.InvalidObsParams):
        val.validate_duration(1e12)


@pytest.mark.parametrize('bad', ['24', True, [1], {'a': 1}, math.nan, math.inf])
def test_duration_bad_types(bad):
    with pytest.raises(val.InvalidObsParams):
        val.validate_duration(bad)


def test_duration_matches_ui_input_max():
    from vlbiplanobs.gui import inputs
    rendered = str(inputs.duration())
    assert f"max={val.MAX_DURATION_H}" in rendered


# ------------------------------------------------------------------------------------------
# Choices from freqsetups
# ------------------------------------------------------------------------------------------
def test_choices_accept_allowed_and_defaults():
    assert val.validate_datarate(None) == val.DEFAULT_DATARATE
    assert val.validate_datarate('1024') == 1024
    for rate in fs.data_rates:
        assert val.validate_datarate(rate) == rate
    assert val.validate_inttime(0.5) == 0.5
    assert val.validate_polarizations(4) == 4


@pytest.mark.parametrize('func,bad', [(val.validate_datarate, 3), (val.validate_datarate, 10**9),
                                      (val.validate_subbands, 0), (val.validate_subbands, 10**6),
                                      (val.validate_channels, 10**8), (val.validate_polarizations, 3),
                                      (val.validate_inttime, 1e-9), (val.validate_inttime, 'abc'),
                                      (val.validate_datarate, [2048]), (val.validate_channels, True)])
def test_choices_reject_unknown(func, bad):
    with pytest.raises(val.InvalidObsParams):
        func(bad)


def test_band_index_and_name():
    assert val.band_name_from_index(0) is None
    assert val.band_name_from_index(None) is None
    assert val.band_name_from_index(1) == list(fs.bands)[0]
    assert val.band_name_from_index(len(fs.bands)) == list(fs.bands)[-1]
    for bad in (len(fs.bands) + 1, -1, '3', 2.0, [1], True):
        with pytest.raises(val.InvalidObsParams):
            val.band_name_from_index(bad)
    with pytest.raises(val.InvalidObsParams):
        val.validate_band('1000cm')


def test_ontarget_percent():
    assert val.ontarget_from_percent(None) == val.DEFAULT_ONTARGET
    assert val.ontarget_from_percent(50) == 0.5
    for bad in (5, 150, 'x', math.nan):
        with pytest.raises(val.InvalidObsParams):
            val.ontarget_from_percent(bad)


def test_start_time_parsing():
    t = val.parse_start_time('2026-10-02T00:00:00', '12:30')
    assert t.datetime.hour == 12 and t.datetime.minute == 30
    for date, time in (('2026-13-40', '00:00'), ('garbage', '00:00'), ('2026-01-01', '25:99'),
                       ('3000-01-01', '00:00'), (None, '00:00'), ('2026-01-01', 5)):
        with pytest.raises(val.InvalidObsParams):
            val.parse_start_time(date, time)


# ------------------------------------------------------------------------------------------
# Targets
# ------------------------------------------------------------------------------------------
def test_clean_targets_filters_and_caps():
    raw = [' 3C273 ', '', '   ', None, 5, {'a': 1}, 'x' * (val.MAX_TARGET_LEN + 1), '3C273', 'M87']
    assert val.clean_target_specs(raw) == ['3C273', 'M87']
    many = [f"src{i}" for i in range(val.MAX_TARGETS * 50)]
    assert val.clean_target_specs(many) == many[:val.MAX_TARGETS]
    assert val.clean_target_specs('3C273') == []
    assert val.clean_target_specs(None) == []


def test_clean_station_codenames():
    assert val.clean_station_codenames(['Ef', 'Ef', 'NotAStation', 3, None]) == ['Ef']
    assert val.clean_station_codenames('Ef') == []


# ------------------------------------------------------------------------------------------
# Whole store-obs-params dict
# ------------------------------------------------------------------------------------------
def test_normalize_obs_params_valid():
    params = val.normalize_obs_params(_valid_params(startdate='2026-10-02', starttime='08:00'))
    assert params['band'] == '18cm' and params['stations'] == ['Ef', 'Mc', 'Nt']
    assert params['start_time'] is not None and params['duration'] == 10.0


@pytest.mark.parametrize('override', [{'duration': 1e9}, {'band': 'xx'}, {'stations': ['Foo']},
                                      {'datarate': 7}, {'inttime': 1e-6}, {'ontarget': 5},
                                      {'startdate': '2026-99-99', 'starttime': '00:00'}])
def test_normalize_obs_params_invalid(override):
    with pytest.raises(val.InvalidObsParams):
        val.normalize_obs_params(_valid_params(**override))


def test_normalize_obs_params_not_a_dict():
    with pytest.raises(val.InvalidObsParams):
        val.normalize_obs_params(['band', '18cm'])


def test_normalize_obs_params_caps_targets():
    params = val.normalize_obs_params(_valid_params(targets=[f"s{i}" for i in range(1000)] + [1, None]))
    assert len(params['targets']) == val.MAX_TARGETS


# ------------------------------------------------------------------------------------------
# url_open-style payload entries
# ------------------------------------------------------------------------------------------
def test_url_component_values():
    group_name, group_stations = next(iter(_groups().items()))
    assert val.validate_url_component({'type': 'group-active-codename', 'index': group_name}, 'data',
                                      sorted(group_stations)[0])[0]
    assert not val.validate_url_component({'type': 'group-active-codename', 'index': group_name}, 'data', 'Ef')[0]
    assert not val.validate_url_component({'type': 'group-active-codename', 'index': group_name}, 'data', 'Nope')[0]
    assert val.validate_url_component('duration', 'value', 12) == (True, 12)
    assert not val.validate_url_component('duration', 'value', 1e9)[0]
    assert not val.validate_url_component('band-slider', 'value', 999)[0]
    assert not val.validate_url_component('datarate', 'value', 123)[0]
    assert not val.validate_url_component('switch-specify-epoch', 'value', 'yes')[0]
    assert not val.validate_url_component('switches-antennas', 'value', ['Ef', 'Nope'])[0]
    assert not val.validate_url_component('switches-antennas', 'value', [['Ef']])[0]
    assert val.validate_url_component('switches-antennas', 'value', ['Ef', 'Mc']) == (True, ['Ef', 'Mc'])
    ok, targets = val.validate_url_component('store-targets', 'data', [f"t{i}" for i in range(500)] + [7])
    assert ok and len(targets) == val.MAX_TARGETS
    assert not val.validate_url_component('starttime', 'value', '99:99')[0]
    assert val.validate_url_component('startdate', 'date', '2026-10-02')[0]
    assert not val.validate_url_component('unknown-component', 'value', 1)[0]
    assert not val.validate_url_component({'type': 'unknown', 'index': 'x'}, 'value', 1)[0]


def _groups() -> dict[str, list[str]]:
    """Return {group: [codenames]} from the station catalogue (skips the test if there are none)."""
    groups: dict[str, list[str]] = {}
    for station in observation._STATIONS:
        if station.group:
            groups.setdefault(station.group, []).append(station.codename)
    if not groups:
        pytest.skip("No grouped stations in the catalogue.")
    return groups


# ------------------------------------------------------------------------------------------
# UV-plot antenna highlight
# ------------------------------------------------------------------------------------------
def test_baseline_has_antenna_exact_match():
    assert plots.baseline_has_antenna('Me-Ef', 'Me')
    assert plots.baseline_has_antenna('Ef-Me', 'Me')
    assert not plots.baseline_has_antenna('Me1-Ef', 'Me')
    assert not plots.baseline_has_antenna('Wb14-Ef', 'Wb')
    assert not plots.baseline_has_antenna('Ef-Mc', 'E')


def test_uvplot_highlight_does_not_match_prefix():
    bl_uv = {'Wb-Ef': np.array([[1.0, 1.0]]), 'Wb14-Ef': np.array([[2.0, 2.0]])}
    fig = plots.uvplot_from_baselines(bl_uv, ['Wb'])
    highlighted = [trace.name for trace in fig.data if trace.marker.color != plots.UV_BASE_COLOR]
    assert highlighted == ['Wb-Ef']


def test_uvplot_one_trace_per_baseline_without_point_text():
    bl_uv = {'Ef-Wb': np.array([[1.0, 2.0], [3.0, 4.0]]), 'Ef-Mc': np.array([[5.0, 6.0]]),
             'Mc-Wb': np.array([[7.0, 8.0]])}
    fig = plots.uvplot_from_baselines(bl_uv, ['Mc'])
    assert sorted(trace.name for trace in fig.data) == sorted(bl_uv)
    # Highlighted traces are drawn last (on top).
    assert [trace.name for trace in fig.data] == ['Ef-Wb', 'Ef-Mc', 'Mc-Wb']
    for trace in fig.data:
        assert trace.text is None and trace.hovertext is None
        assert '%{fullData.name}' in trace.hovertemplate
    efwb = fig.data[0]
    assert list(efwb.x) == [1.0, 3.0, -1.0, -3.0] and list(efwb.y) == [2.0, 4.0, -2.0, -4.0]
    assert efwb.marker.color == plots.UV_BASE_COLOR and efwb.marker.size == plots.UV_BASE_SIZE
    assert fig.data[1].marker.color == plots.UV_HIGHLIGHT_COLORS[0]
    assert fig.data[1].marker.size == plots.UV_HIGHLIGHT_SIZE
    # Serialized figure carries the uv arrays as compact base64 typed arrays, not per-point lists.
    payload = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    assert '"bdata"' in payload and 'Ef-Wb' in payload


def test_uv_marker_style_color_follows_selection_order():
    assert plots.uv_marker_style('Ef-Wb', None)['highlighted'] is False
    assert plots.uv_marker_style('Ef-Wb', ['Mc', 'Wb'])['color'] == plots.UV_HIGHLIGHT_COLORS[1]
    assert plots.uv_marker_style('Ef-Wb', ['Ef', 'Wb'])['color'] == plots.UV_HIGHLIGHT_COLORS[0]


def test_uv_highlight_is_clientside_and_store_removed():
    from vlbiplanobs.gui import callbacks
    import inspect
    from vlbiplanobs.gui import outputs
    assert 'store-uv-data' not in inspect.getsource(outputs)
    assert 'store-uv-data' not in inspect.getsource(callbacks)
    assert "split('-')" in callbacks.uv_highlight_javascript


# ------------------------------------------------------------------------------------------
# Polaris origin allowlist (PLANOBS_POLARIS_ORIGINS)
# ------------------------------------------------------------------------------------------
@pytest.mark.parametrize('raw,expected', [(None, []), ('', []), (' , ', []),
                                          ('https://polaris.example.org', ['https://polaris.example.org']),
                                          ('HTTPS://Polaris.Example.org:443/', ['https://polaris.example.org']),
                                          ('http://localhost:8080, http://localhost:80',
                                           ['http://localhost:8080', 'http://localhost']),
                                          ('https://a.org,https://a.org', ['https://a.org']),
                                          ('https://[::1]:8443', ['https://[::1]:8443'])])
def test_parse_polaris_origins(raw, expected):
    assert val.parse_polaris_origins(raw) == expected


@pytest.mark.parametrize('bad', ['*', 'polaris.example.org', 'ftp://a.org', 'https://a.org/path', 'https://a.org?x=1',
                                 'https://user@a.org', 'https://a.org:99999', 'javascript:alert(1)', 'null'])
def test_parse_polaris_origins_rejects_invalid(bad):
    assert val.parse_polaris_origins(f"{bad},https://ok.org") == ['https://ok.org']


def test_polaris_export_never_posts_to_wildcard():
    from vlbiplanobs.gui import callbacks
    assert "postMessage(value, '*')" not in callbacks.callback_javascript
    assert 'postMessage(value, polarisOrigin)' in callbacks.callback_javascript
    assert 'document.referrer' in callbacks.callback_javascript
