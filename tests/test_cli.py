import pytest
from astropy import units as u
from vlbiplanobs.cli import VLBIObs
from vlbiplanobs import observation as obs

def test_vlbiobs_init_minimal():
    # Minimal valid VLBIObs initialization
    o = VLBIObs(band='18cm', stations=obs._STATIONS.filter_antennas(['Ef']), scans={})
    assert o.band == '18cm'
    assert hasattr(o, 'stations')
    assert hasattr(o, 'scans')

def test_vlbiobs_inherits_observation():
    # VLBIObs should be a subclass of Observation
    o = VLBIObs(band='18cm', stations=obs._STATIONS.filter_antennas(['Ef']), scans={})
    assert isinstance(o, obs.Observation)

def test_vlbiobs_summary_methods():
    o = VLBIObs(band='18cm', stations=obs._STATIONS.filter_antennas(['Ef']), scans={})
    # These methods should run without error
    o.summary(gui=False, tui=True)
    o.summary(gui=True, tui=False)

def test_vlbiobs_properties():
    evn6 = ['Ef', 'Jb2', 'O8', 'T6', 'Wb']
    o = VLBIObs(band='18cm', stations=obs._STATIONS.filter_antennas(evn6), duration=10*u.h, scans={})
    assert o.duration == 10*u.h
    assert set(o.stations.station_codenames) == set(evn6)

def test_vlbiobs_invalid_init():
    # Should raise if required arguments are missing
    with pytest.raises(TypeError):
        VLBIObs()
