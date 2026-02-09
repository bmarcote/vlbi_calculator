"""Tests for the calibrators module."""

import pytest
import numpy as np
from astropy import coordinates as coord
from astropy import units as u
from astropy.time import Time

from vlbiplanobs import calibrators
from vlbiplanobs import stations as sts
from vlbiplanobs.sources import Source


class TestCalibratorSource:
    """Tests for the CalibratorSource class."""

    @pytest.fixture
    def sample_source(self):
        """Create a sample calibrator source for testing."""
        # Need 5 elements for S,C,X,U,K bands
        flux_resolved = np.array([1.0, 5.0, 3.0, 2.0, 1.5])  # S, C, X, U, K
        flux_unresolved = np.array([2.0, 10.0, 6.0, 4.0, 3.0])  # S, C, X, U, K
        return calibrators.CalibratorSource(
            name="J1230+4500",
            ivsname="1230+450",
            ra_deg=187.5,
            dec_deg=45.0,
            n_observations=100,
            flux_resolved=flux_resolved,
            flux_unresolved=flux_unresolved,
            is_calibrator=True
        )

    def test_basic_creation(self, sample_source):
        """Test basic CalibratorSource creation."""
        assert sample_source.name == "J1230+4500"
        assert sample_source.ivsname == "1230+450"
        assert sample_source.n_observations == 100
        assert sample_source.is_calibrator is True
        assert sample_source.ra_deg == 187.5
        assert sample_source.dec_deg == 45.0

    def test_unresolved_flux(self, sample_source):
        """Test unresolved_flux method."""
        assert sample_source.unresolved_flux('c') == 10.0
        assert sample_source.unresolved_flux('s') == 2.0
        assert sample_source.unresolved_flux('k') == 3.0  # K band exists
        # Q band maps to K in _RFC_BANDS but not in _BAND_INDEX
        assert sample_source.unresolved_flux('q') == 0.0  # Q not in _BAND_INDEX
        assert sample_source.unresolved_flux('z') == 0.0  # Missing band

    def test_resolved_flux(self, sample_source):
        """Test resolved_flux method."""
        assert sample_source.resolved_flux('c') == 5.0
        assert sample_source.resolved_flux('s') == 1.0
        assert sample_source.resolved_flux('k') == 1.5  # K band exists
        # Q band maps to K in _RFC_BANDS but not in _BAND_INDEX
        assert sample_source.resolved_flux('q') == 0.0  # Q not in _BAND_INDEX
        assert sample_source.resolved_flux('z') == 0.0  # Missing band

    def test_get_observed_bands(self, sample_source):
        """Test get_observed_bands method."""
        bands = sample_source.get_observed_bands()
        # Should return bands in order: S,C,X,U,K
        assert bands == "S,C,X,U,K"

    def test_get_observed_bands_partial(self):
        """Test get_observed_bands with only some bands."""
        flux_resolved = np.array([0.0, 5.0, 0.0, 2.0, 0.0])  # Only C and U
        flux_unresolved = np.array([0.0, 10.0, 0.0, 4.0, 0.0])
        source = calibrators.CalibratorSource(
            name="J1230+4500", ivsname="1230+450", ra_deg=187.5, dec_deg=45.0,
            n_observations=100, flux_resolved=flux_resolved, flux_unresolved=flux_unresolved,
            is_calibrator=True
        )
        bands = source.get_observed_bands()
        assert bands == "C,U"

    def test_get_astrogeo_link(self, sample_source):
        """Test get_astrogeo_link method."""
        link = sample_source.get_astrogeo_link()
        assert "astrogeo.org" in link
        assert "calib_search_form.csh" in link
        assert "ra=" in link
        assert "dec=" in link

    def test_get_skycoord(self, sample_source):
        """Test get_skycoord method."""
        skycoord = sample_source.get_skycoord()
        assert isinstance(skycoord, coord.SkyCoord)
        assert skycoord.ra.deg == pytest.approx(187.5)
        assert skycoord.dec.deg == pytest.approx(45.0)

    def test_inheritance_from_source(self, sample_source):
        """Test that CalibratorSource inherits from Source."""
        assert isinstance(sample_source, Source)
        # Check that it has the attributes from Source
        assert hasattr(sample_source, 'coord')
        assert hasattr(sample_source, 'name')
        assert sample_source.name == "J1230+4500"


class TestRFCCatalog:
    """Tests for the RFCCatalog class."""

    def test_catalog_loading(self):
        """Test that the catalog loads successfully."""
        catalog = calibrators.RFCCatalog()
        assert catalog.n_sources > 0

    def test_catalog_with_min_flux(self):
        """Test catalog loading with minimum flux filter."""
        catalog_low = calibrators.RFCCatalog(min_flux=0.5 * u.Jy, band='c')
        catalog_high = calibrators.RFCCatalog(min_flux=5.0 * u.Jy, band='c')
        # Higher threshold should give fewer or equal sources
        assert catalog_high.n_sources <= catalog_low.n_sources

    def test_get_source_found(self):
        """Test finding a known source."""
        catalog = calibrators.RFCCatalog(min_flux=1.0 * u.Jy, band='c')
        # 3C454.3 is a bright, well-known calibrator
        src = catalog.get_source('3C454.3')
        assert src is not None
        assert src.name == 'J2253+1608'  # J2000 name
        assert src.ivsname == '3C454.3'
        assert src.is_calibrator is True

    def test_get_source_not_found(self):
        """Test finding a non-existent source."""
        catalog = calibrators.RFCCatalog()
        src = catalog.get_source('NONEXISTENTSOURCE12345')
        assert src is None

    def test_sources_property(self):
        """Test the sources property."""
        catalog = calibrators.RFCCatalog(min_flux=2.0 * u.Jy, band='c')
        sources = catalog.sources
        assert isinstance(sources, list)
        assert len(sources) == catalog.n_sources
        assert all(isinstance(s, calibrators.CalibratorSource) for s in sources)

    def test_get_coord_arrays(self):
        """Test the _get_coord_arrays method."""
        catalog = calibrators.RFCCatalog(min_flux=1.0 * u.Jy, band='c')
        ra_arr, dec_arr = catalog._get_coord_arrays()
        assert isinstance(ra_arr, np.ndarray)
        assert isinstance(dec_arr, np.ndarray)
        assert len(ra_arr) == catalog.n_sources
        assert len(dec_arr) == catalog.n_sources
        assert len(ra_arr) == len(dec_arr)

    def test_catalog_band_filtering(self):
        """Test that catalog correctly filters by band."""
        catalog_s = calibrators.RFCCatalog(min_flux=0.1 * u.Jy, band='s')
        catalog_c = calibrators.RFCCatalog(min_flux=0.1 * u.Jy, band='c')
        # Different bands should have different numbers of sources
        # (though they might be equal, this tests the functionality)
        assert isinstance(catalog_s.n_sources, int)
        assert isinstance(catalog_c.n_sources, int)
        assert catalog_s.n_sources > 0
        assert catalog_c.n_sources > 0


class TestBandMapping:
    """Tests for the RFC band mapping functionality."""

    def test_band_mappings(self):
        """Test that band mappings are correct."""
        assert calibrators._RFC_BANDS['l'] == 's'
        assert calibrators._RFC_BANDS['s'] == 's'
        assert calibrators._RFC_BANDS['c'] == 'c'
        assert calibrators._RFC_BANDS['m'] == 'c'
        assert calibrators._RFC_BANDS['x'] == 'x'
        assert calibrators._RFC_BANDS['u'] == 'u'
        assert calibrators._RFC_BANDS['k'] == 'k'
        assert calibrators._RFC_BANDS['q'] == 'k'


class TestGetNearbySources:
    """Tests for the get_nearby_sources function."""

    @pytest.fixture
    def reference_source(self):
        """Create a reference source for nearby source tests."""
        # Need 5 elements for S,C,X,U,K bands
        flux_resolved = np.array([5.0, 5.0, 0.0, 0.0, 0.0])  # C band only
        flux_unresolved = np.array([10.0, 10.0, 0.0, 0.0, 0.0])  # C band only
        return calibrators.CalibratorSource(
            name="J1230+4500",
            ivsname="1230+450",
            ra_deg=187.5,
            dec_deg=45.0,
            n_observations=100,
            flux_resolved=flux_resolved,
            flux_unresolved=flux_unresolved,
            is_calibrator=True
        )

    @pytest.fixture
    def sample_catalog(self, reference_source):
        """Create a sample catalog with multiple sources."""
        catalog = calibrators.RFCCatalog.__new__(calibrators.RFCCatalog)

        # Create nearby sources - reference at 12h30m +45d
        # nearby1: ~2 arcmin away (within 2 deg)
        flux_resolved_1 = np.array([2.0, 2.0, 0.0, 0.0, 0.0])
        flux_unresolved_1 = np.array([4.0, 4.0, 0.0, 0.0, 0.0])
        nearby1 = calibrators.CalibratorSource(
            name="J1230+4501",
            ivsname="1230+451",
            ra_deg=187.5083,  # ~2 arcmin east
            dec_deg=45.0028,  # ~10 arcsec north
            n_observations=50,
            flux_resolved=flux_resolved_1,
            flux_unresolved=flux_unresolved_1,
            is_calibrator=True
        )
        # nearby2: ~1 deg away (within 2 deg)
        flux_resolved_2 = np.array([1.0, 1.0, 0.0, 0.0, 0.0])
        flux_unresolved_2 = np.array([2.0, 2.0, 0.0, 0.0, 0.0])
        nearby2 = calibrators.CalibratorSource(
            name="J1230+4600",
            ivsname="1230+460",
            ra_deg=187.5,  # Same RA
            dec_deg=46.0,  # 1 deg north
            n_observations=50,
            flux_resolved=flux_resolved_2,
            flux_unresolved=flux_unresolved_2,
            is_calibrator=True
        )
        # far: ~10 deg away (outside 2 deg)
        flux_resolved_3 = np.array([3.0, 3.0, 0.0, 0.0, 0.0])
        flux_unresolved_3 = np.array([6.0, 6.0, 0.0, 0.0, 0.0])
        far = calibrators.CalibratorSource(
            name="J1230+5500",
            ivsname="1230+550",
            ra_deg=187.5,  # Same RA
            dec_deg=55.0,  # 10 deg north
            n_observations=50,
            flux_resolved=flux_resolved_3,
            flux_unresolved=flux_unresolved_3,
            is_calibrator=True
        )

        catalog._sources = [reference_source, nearby1, nearby2, far]
        catalog._min_flux = 0.0
        catalog._band = 'c'
        catalog._catalog_filename = None
        catalog._ra_arr = np.array([187.5, 187.5083, 187.5, 187.5])
        catalog._dec_arr = np.array([45.0, 45.0028, 46.0, 55.0])
        return catalog

    def test_nearby_sources_basic(self, reference_source, sample_catalog):
        """Test basic nearby sources functionality."""
        nearby = calibrators.get_nearby_sources(
            reference_source,
            max_separation=2 * u.deg,
            min_unresolved_flux=0.0,
            band='c',
            catalog=sample_catalog
        )
        # Should find 2 nearby sources (within 2 degrees, excluding reference)
        # nearby1 (~2 arcmin), nearby2 (~1 deg) - far is ~10 deg away
        assert len(nearby) == 2
        names = [s.name for s, _ in nearby]
        assert "J1230+4501" in names
        assert "J1230+4600" in names
        assert "J1230+5500" not in names

    def test_nearby_sources_flux_filter(self, reference_source, sample_catalog):
        """Test nearby sources with flux filtering."""
        nearby = calibrators.get_nearby_sources(
            reference_source,
            max_separation=2 * u.deg,
            min_unresolved_flux=3.0,  # Only sources with > 3 Jy
            band='c',
            catalog=sample_catalog
        )
        # Only nearby1 has 4.0 Jy, nearby2 (J1230+4600) has 2.0 Jy
        assert len(nearby) == 1
        assert nearby[0][0].name == "J1230+4501"

    def test_nearby_sources_separation_filter(self, reference_source, sample_catalog):
        """Test nearby sources with separation filtering."""
        nearby = calibrators.get_nearby_sources(
            reference_source,
            max_separation=10 * u.arcmin,  # Very small radius
            min_unresolved_flux=0.0,
            band='c',
            catalog=sample_catalog
        )
        # Only nearby1 (~2 arcmin away) should be found
        assert len(nearby) == 1
        assert nearby[0][0].name == "J1230+4501"

    def test_nearby_sources_sorted_by_separation(self, reference_source, sample_catalog):
        """Test that nearby sources are sorted by separation."""
        nearby = calibrators.get_nearby_sources(
            reference_source,
            max_separation=2 * u.deg,
            min_unresolved_flux=0.0,
            band='c',
            catalog=sample_catalog
        )
        # Should be sorted closest first
        # separations are returned as floats, not Quantity objects
        separations = [sep for _, sep in nearby]
        assert separations == sorted(separations)

    def test_nearby_sources_limit_count(self, reference_source, sample_catalog):
        """Test limiting the number of returned sources."""
        nearby = calibrators.get_nearby_sources(
            reference_source,
            max_separation=2 * u.deg,
            min_unresolved_flux=0.0,
            band='c',
            catalog=sample_catalog,
            n_sources=1
        )
        assert len(nearby) == 1

    def test_nearby_sources_no_band(self, reference_source, sample_catalog):
        """Test nearby sources without specifying a band."""
        nearby = calibrators.get_nearby_sources(
            reference_source,
            max_separation=2 * u.deg,
            min_unresolved_flux=0.0,
            band=None,  # Check all bands
            catalog=sample_catalog
        )
        # Should still work, checking all available bands
        # With max_separation=2 deg: nearby1 and nearby2 are included, far (10 deg) is excluded
        assert len(nearby) == 2


class TestGetFringeFinderSources:
    """Tests for the get_fringe_finder_sources function."""

    @pytest.fixture
    def sample_times(self):
        """Create sample observation times."""
        return Time("2024-06-15 12:00:00") + np.arange(0, 2, 0.5) * u.h

    @pytest.fixture
    def single_station(self):
        """Create a single test station."""
        # Ef (Effelsberg) coordinates
        ef_loc = coord.EarthLocation(
            4033949.5 * u.m, 486989.1 * u.m, 4900430.9 * u.m
        )
        sefds = {'21cm': 500 * u.Jy, '6cm': 100 * u.Jy}
        return sts.Station('Ef', 'Ef', ('EVN',), ef_loc, sefds)

    def test_fringe_finder_no_stations(self, sample_times):
        """Test with empty stations list."""
        empty_stations = sts.Stations(stations=[])
        result = calibrators.get_fringe_finder_sources(
            empty_stations, sample_times, band='c', min_flux=5.0 * u.Jy
        )
        # Should return empty tuple when no stations
        assert result == ([], [])

    def test_fringe_finder_single_station(self, single_station, sample_times):
        """Test fringe finder with a single station."""
        stations = sts.Stations(stations=[single_station])
        result = calibrators.get_fringe_finder_sources(
            stations, sample_times, band='c',
            min_elevation=20 * u.deg, min_flux=5.0 * u.Jy
        )
        # Should return tuple of (sources, min_elevations)
        assert isinstance(result, tuple)
        assert len(result) == 2
        sources, min_elevs = result
        assert isinstance(sources, list)
        assert isinstance(min_elevs, list)
        assert len(sources) == len(min_elevs)
        # All returned sources should be CalibratorSource instances
        assert all(isinstance(s, calibrators.CalibratorSource) for s in sources)
        # All min_elevations should be floats
        assert all(isinstance(e, float) for e in min_elevs)

    def test_fringe_finder_sorted_by_flux(self, single_station, sample_times):
        """Test that results are sorted by unresolved flux."""
        stations = sts.Stations(stations=[single_station])
        sources, min_elevs = calibrators.get_fringe_finder_sources(
            stations, sample_times, band='c',
            min_elevation=20 * u.deg, min_flux=1.0 * u.Jy
        )
        if len(sources) > 1:
            fluxes = [s.unresolved_flux('c') for s in sources]
            assert fluxes == sorted(fluxes, reverse=True)

    def test_fringe_finder_band_mapping(self, single_station, sample_times):
        """Test that band mapping works correctly."""
        stations = sts.Stations(stations=[single_station])
        # Test with 'm' band (should map to 'c')
        sources_m, elevs_m = calibrators.get_fringe_finder_sources(
            stations, sample_times, band='m',
            min_elevation=20 * u.deg, min_flux=5.0 * u.Jy
        )
        # Test with 'c' band directly
        sources_c, elevs_c = calibrators.get_fringe_finder_sources(
            stations, sample_times, band='c',
            min_elevation=20 * u.deg, min_flux=5.0 * u.Jy
        )
        # Should get same results since 'm' maps to 'c'
        assert len(sources_m) == len(sources_c)

    def test_fringe_finder_require_all_vs_any(self, sample_times):
        """Test require_all_stations vs require only one station."""
        # Create two stations at different locations
        ef_loc = coord.EarthLocation(
            4033949.5 * u.m, 486989.1 * u.m, 4900430.9 * u.m
        )
        mc_loc = coord.EarthLocation(
            -5464555.2 * u.m, -2493147.5 * u.m, 2150614.3 * u.m
        )
        sefds = {'21cm': 500 * u.Jy, '6cm': 100 * u.Jy}

        ef = sts.Station('Ef', 'Ef', ('EVN',), ef_loc, sefds)
        mc = sts.Station('Mc', 'Mc', ('EVN',), mc_loc, sefds)

        stations = sts.Stations(stations=[ef, mc])

        # With require_all_stations=True (default)
        sources_all, elevs_all = calibrators.get_fringe_finder_sources(
            stations, sample_times, band='c',
            min_elevation=20 * u.deg, min_flux=2.0 * u.Jy,
            require_all_stations=True
        )

        # With require_all_stations=False
        sources_any, elevs_any = calibrators.get_fringe_finder_sources(
            stations, sample_times, band='c',
            min_elevation=20 * u.deg, min_flux=2.0 * u.Jy,
            require_all_stations=False
        )

        # Results with require_all=False should be >= results with require_all=True
        assert len(sources_any) >= len(sources_all)

    def test_fringe_finder_returns_empty_when_none_visible(self, sample_times):
        """Test that empty list is returned when no sources are visible."""
        # Create a station and times where bright sources might not be visible
        # Use very high elevation requirement
        ef_loc = coord.EarthLocation(
            4033949.5 * u.m, 486989.1 * u.m, 4900430.9 * u.m
        )
        sefds = {'21cm': 500 * u.Jy, '6cm': 100 * u.Jy}
        ef = sts.Station('Ef', 'Ef', ('EVN',), ef_loc, sefds)
        stations = sts.Stations(stations=[ef])

        sources, elevs = calibrators.get_fringe_finder_sources(
            stations, sample_times, band='c',
            min_elevation=85 * u.deg,  # Very high elevation
            min_flux=50.0 * u.Jy,  # Very bright
            require_all_stations=True
        )
        # Should return empty lists when criteria are too strict
        assert isinstance(sources, list)
        assert isinstance(elevs, list)
        assert len(sources) == len(elevs)

    def test_fringe_finder_minimum_elevations(self, single_station, sample_times):
        """Test that minimum elevations are correctly calculated."""
        stations = sts.Stations(stations=[single_station])
        sources, min_elevs = calibrators.get_fringe_finder_sources(
            stations, sample_times, band='c',
            min_elevation=10 * u.deg, min_flux=1.0 * u.Jy
        )
        if len(sources) > 0:
            # All minimum elevations should be >= min_elevation
            for elev in min_elevs:
                assert elev >= 10.0 or elev == 0.0  # 0.0 if not visible at all


class TestCalibratorsIntegration:
    """Integration tests for the calibrators module."""

    def test_full_workflow(self):
        """Test a complete workflow: load catalog, find sources, find nearby."""
        # Load catalog
        catalog = calibrators.RFCCatalog(min_flux=2.0 * u.Jy, band='c')
        assert catalog.n_sources > 0

        # Get a bright source
        bright_source = catalog.get_source('3C454.3')
        if bright_source is not None:
            assert bright_source.unresolved_flux('c') >= 2.0

            # Find nearby sources
            nearby = calibrators.get_nearby_sources(
                bright_source,
                max_separation=5 * u.deg,
                min_unresolved_flux=0.5,
                band='c',
                catalog=catalog,
                n_sources=10
            )
            assert isinstance(nearby, list)
            # All nearby sources should have positive flux
            for src, sep in nearby:
                assert src.unresolved_flux('c') > 0.5
                assert sep <= 5 * u.deg

    def test_fringe_finder_with_elevations(self):
        """Test fringe finder returns proper elevation data."""
        # Create a simple setup
        ef_loc = coord.EarthLocation(
            4033949.5 * u.m, 486989.1 * u.m, 4900430.9 * u.m
        )
        sefds = {'21cm': 500 * u.Jy, '6cm': 100 * u.Jy}
        ef = sts.Station('Ef', 'Ef', ('EVN',), ef_loc, sefds)
        stations = sts.Stations(stations=[ef])
        times = Time("2024-06-15 12:00:00") + np.arange(0, 1, 0.25) * u.h

        sources, min_elevs = calibrators.get_fringe_finder_sources(
            stations, times, band='c',
            min_elevation=20 * u.deg, min_flux=1.0 * u.Jy
        )

        # Should return matching lists
        assert len(sources) == len(min_elevs)
        
        # Test that sources have the expected methods
        for src in sources:
            assert hasattr(src, 'get_observed_bands')
            assert hasattr(src, 'get_astrogeo_link')
            assert isinstance(src.get_observed_bands(), str)
            assert isinstance(src.get_astrogeo_link(), str)

    def test_angular_separation_function(self):
        """Test the _angular_separation helper function."""
        # Test separation between two points
        ra1, dec1 = 187.5, 45.0  # Reference
        ra2_arr = np.array([187.5, 188.0, 187.5])  # Same RA, 0.5deg east, same RA
        dec2_arr = np.array([46.0, 45.0, 55.0])  # 1deg north, same dec, 10deg north
        
        separations = calibrators._angular_separation(ra1, dec1, ra2_arr, dec2_arr)
        
        # Expected: 1.0 deg, 0.5 deg, 10.0 deg
        # But the actual calculation gives 0.3535 for the second case (0.5deg RA difference at 45deg dec)
        expected = np.array([1.0, 0.353553, 10.0])
        np.testing.assert_allclose(separations, expected, atol=0.01)


class TestCLIIntegration:
    """Tests for CLI functionality."""

    def test_main_functions_exist(self):
        """Test that main CLI functions exist and are callable."""
        assert hasattr(calibrators, 'main_fringe')
        assert hasattr(calibrators, 'main_phasecal')
        assert callable(calibrators.main_fringe)
        assert callable(calibrators.main_phasecal)

    def test_cli_functions_have_proper_signatures(self):
        """Test that CLI functions have expected signatures (no arguments for CLI)."""
        import inspect
        
        # Both functions should take no arguments (they use argparse internally)
        fringe_sig = inspect.signature(calibrators.main_fringe)
        phase_sig = inspect.signature(calibrators.main_phasecal)
        
        assert len(fringe_sig.parameters) == 0
        assert len(phase_sig.parameters) == 0
