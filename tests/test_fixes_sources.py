"""Regression tests for fixes in sources.py, calibrators.py and constraints.py."""
import numpy as np
import pytest
from astropy import units as u
from astropy import coordinates as coord
from vlbiplanobs import sources, calibrators


def _src(name: str, src_type: sources.SourceType) -> sources.Source:
    """Returns a Source at fixed coordinates with the given name and type."""
    return sources.Source(name, '10h20m10s 40d30m10s', source_type=src_type)


# --- 1. RFCCatalog.get_source exact matching --------------------------------------------------------

def test_get_source_no_substring_match():
    catalog = calibrators.RFCCatalog(min_flux=0.0 * u.Jy, band='c')
    assert catalog.get_source('1') is None
    assert catalog.get_source('J0000') is None
    found = catalog.get_source('j2253+1608')
    assert found is not None and found.name == 'J2253+1608'
    assert catalog.get_source('3c454.3').name == 'J2253+1608'


def test_get_source_filtered_out_returns_none():
    full = calibrators.RFCCatalog(min_flux=0.0 * u.Jy, band='c')
    bright = calibrators.RFCCatalog(min_flux=5.0 * u.Jy, band='c')
    bright_names = {s.name for s in bright.sources}
    faint = next(s for s in full.sources if s.name not in bright_names)
    assert bright.get_source(faint.name) is None
    assert bright.get_source(faint.ivsname) is None


# --- 2. Source kwargs / SkyCoord reuse / catalog isolation -------------------------------------------

def test_source_kwargs_passed_to_skycoord():
    src = sources.Source('x', '187.7 12.39', unit='deg')
    assert src.coord.ra.deg == pytest.approx(187.7)
    assert src.coord.dec.deg == pytest.approx(12.39)


def test_source_reuses_skycoord():
    sky = coord.SkyCoord(ra=10 * u.deg, dec=20 * u.deg)
    assert sources.Source('x', sky).coord is sky


def test_rfc_catalog_instances_independent():
    cat1 = calibrators.RFCCatalog(min_flux=1.0 * u.Jy, band='c')
    cat2 = calibrators.RFCCatalog(min_flux=1.0 * u.Jy, band='c')
    assert cat1.sources is not cat2.sources
    assert cat1.n_sources == cat2.n_sources
    cat1.sources.pop()
    assert cat2.n_sources == cat1.n_sources + 1
    with pytest.raises(ValueError):
        cat2.sources[0].flux_unresolved[0] = 99.0
    src = cat2.sources[0]
    assert src.coord.ra.deg == pytest.approx(src.ra_deg)
    assert np.allclose(cat2._get_coord_arrays()[0], [s.ra_deg for s in cat2.sources])


def test_rfc_catalog_invalid_band():
    with pytest.raises(ValueError):
        calibrators.RFCCatalog(band='z')


# --- 4. Name validation before online resolution -------------------------------------------------------

@pytest.mark.parametrize('bad_name', ['', 'a;b', 'x' * 81, 'name<script>', 'a&b=c'])
def test_resolve_name_online_rejects_invalid(bad_name):
    with pytest.raises(ValueError):
        sources.resolve_name_online(bad_name)


def test_validate_source_name_accepts_common_names():
    for name in ('Cyg X-1', 'PSR J0437-4715', "Barnard's star", 'NGC 1068', '3C 84', 'M 87*'):
        assert sources.validate_source_name(name) == name


# --- 5. astrogeo link ---------------------------------------------------------------------------------

def test_astrogeo_link_small_negative_dec():
    flux = np.zeros(5)
    src = calibrators.CalibratorSource('J0000-0030', '2357-007', 0.0, -0.5, 1, flux, flux, True)
    link = src.get_astrogeo_link()
    assert link.startswith('https://')
    assert 'dec=-00%3A30%3A00.000' in link


def test_astrogeo_link_positive_dec():
    flux = np.zeros(5)
    src = calibrators.CalibratorSource('J0000+0030', '2357+007', 0.0, 0.5, 1, flux, flux, True)
    assert 'dec=%2B00%3A30%3A00.000' in src.get_astrogeo_link()


# --- 6 & 8. TOML catalog parsing ----------------------------------------------------------------------

def test_toml_without_duration_uses_default(tmp_path):
    toml_file = tmp_path / 'cat.toml'
    toml_file.write_text('[[target]]\nname = "T1"\ncoordinates = {RA = "10:20:10", Dec = "40:30:10"}\n'
                         '[target.phasecal]\nname = "P1"\ncoordinates = {RA = "10:21:00", Dec = "40:31:00"}\n')
    catalog = sources.SourceCatalog(str(toml_file))
    block = catalog.targets['T1']
    assert all(s.duration == 10 * u.min for s in block.scans)
    assert len(block.fill(2 * u.h)) > 0


def test_pulsar_only_catalog(tmp_path):
    toml_file = tmp_path / 'psr.toml'
    toml_file.write_text('[[pulsar]]\nname = "PSR1"\nduration = 5\n'
                         'coordinates = {RA = "10:20:10", Dec = "-40:30:10"}\n')
    catalog = sources.SourceCatalog(str(toml_file))
    assert catalog.targets == {}
    assert catalog.source_names() == []
    assert catalog.sources() == {}
    assert catalog.source_names(include_calibrators=True) == ['PSR1']
    assert catalog.pulsars['PSR1'].scans[0].duration == 5 * u.min


# --- 7. Scan / ScanBlock.fill -----------------------------------------------------------------------

def test_scan_every_zero_raises():
    with pytest.raises(ValueError):
        sources.Scan(_src('c', sources.SourceType.CHECKSOURCE), every=0)
    with pytest.raises(ValueError):
        sources.Scan(_src('c', sources.SourceType.CHECKSOURCE), every=-3)


def test_fill_no_phasecal_with_periodic_check():
    target = sources.Scan(_src('T', sources.SourceType.TARGET), duration=2 * u.min)
    check = sources.Scan(_src('C', sources.SourceType.CHECKSOURCE), duration=2 * u.min, every=3)
    names = [s.source.name for s in sources.ScanBlock([target, check]).fill(12 * u.min)]
    assert names == ['T', 'T', 'C', 'T', 'T', 'C']


def test_fill_zero_duration_raises():
    target = sources.Scan(_src('T', sources.SourceType.TARGET), duration=0 * u.min)
    pcal = sources.Scan(_src('P', sources.SourceType.PHASECAL), duration=0 * u.min)
    with pytest.raises(ValueError):
        sources.ScanBlock([pcal, target]).fill(1 * u.h)
    target2 = sources.Scan(_src('T', sources.SourceType.TARGET), duration=1 * u.min)
    check = sources.Scan(_src('C', sources.SourceType.CHECKSOURCE), duration=0 * u.min, every=1)
    with pytest.raises(ValueError):
        sources.ScanBlock([target2, check]).fill(1 * u.h)


def test_fractional_time_reflects_mutation():
    target = sources.Scan(_src('T', sources.SourceType.TARGET), duration=2 * u.min)
    pcal = sources.Scan(_src('P', sources.SourceType.PHASECAL), duration=2 * u.min)
    block = sources.ScanBlock([pcal, target])
    assert float(block.fractional_time()['T']) == pytest.approx(0.5)
    target.duration = 6 * u.min
    assert float(block.fractional_time()['T']) == pytest.approx(0.75)
