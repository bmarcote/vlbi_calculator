import pytest
from functools import reduce
import operator
import numpy as np
from astropy import coordinates as coord
from astropy import units as u
from astropy.time import Time
from rich import print as rprint
from vlbiplanobs import sources
# from vlbiplanobs import stations


def test_source():
    """Tests the Source class
    """
    with pytest.raises(TypeError):
        sources.Source()

    with pytest.raises(ValueError):
        sources.Source('a_name', 'coordinates')
        sources.Source('a_name', '00h00m00s xxhqqmees')
        sources.Source('a_name', '40h00m00s 10d10m10s')
        sources.Source('a_name', '00:00:00 100:00:00')
        sources.Source('a_name', '-30:00:00 10:00:00')
        # sources.Source('Cyg X-1')
        sources.Source(name='WQQWE')
        sources.Source('WQQWE')

    s1 = sources.Source('a_name', '10h20m10s 40d30m10s')
    s2 = sources.Source('a_name', '10:20:10 40:30:10', unit=(u.hourangle, u.deg))

    with pytest.raises(ValueError):
        sources.Source(['a_name', 'a_second_name'])
        sources.Source(['a_name', 'a_second_name'], ['10h20m10s 40d30m10s', '10h20m10s 40d30m10s'])
        sources.Source(coordinates=['10h20m10s 40d30m10s', '20h30m0s 20d20m10s'])

    try:
        s3 = sources.Source(name='Cyg X-1')
        assert s3.name == 'Cyg X-1'
        s3.coord.to_string('hmsdms')
        assert s3.coord.separation(coord.SkyCoord('19h58m21.6757344s +35d12m05.784516s')) < 0.01*u.arcsec
    except coord.name_resolve.NameResolveError:
        rprint("[yellow]WARNING: seems like no internet connection is available. "
               "Unable to resolve source names[/yellow]")

    assert s1.name == 'a_name'
    assert s2.name == 'a_name'
    assert s1.coord == s2.coord

    # checking that the distance to the Sun is right
    # times are spring equinox, summer solstice fall equinox
    times = Time(['2024-03-19 22:06', '2024-06-20 15:51', '2024-09-22 07:44'], scale='utc')
    s1 = sources.Source('equinox source', '0h0m0s 0d0m0s')
    # print(f"\n\nSeparations: {s1.sun_separation(times)}\n\n")
    assert all(np.abs(s1.sun_separation(times) - np.array([0.0, 90.0, 180])*u.deg) < 1*u.deg)
    # print(f"\n\nConstraints: {s1.sun_constraint(20*u.deg, times=times)}\n\n")
    assert len(s1.sun_constraint(20*u.deg, times=times)) == 1
    s1 = sources.Source('polar_source', '0h0m0s 90d0m0s')
    # print(f"\n\nSeparations: {s1.sun_separation(times)}\n\n")
    assert all(np.abs(s1.sun_separation(times[(0, 2),]) - np.array([90.0, 90.0])*u.deg) < 1*u.deg)

    # Playing with multiple sources
    # s1 = sources.Source(['s2', 's3'], ['10h20m10s 40d30m10s', '20h30m0s 20d20m10s'])
    # sep = s1.sun_separation(times)
    # assert len(sep) == 2
    # assert (len(sep[0]) == len(times)) and (len(sep[1]) == len(times))
    # s1.sun_constraint(20*u.deg)


def test_source_from_names():
    try:
        assert isinstance(sources.Source.source_from_name('Cyg X-1'), sources.Source)
        assert isinstance(sources.Source.get_coordinates_from_name('Cyg X-1'), coord.SkyCoord)

        src_name1 = '2358+390'
        src_name2 = 'J0000+3918'
        cor1 = sources.Source.get_coordinates_from_name(src_name1)
        cor2 = sources.Source.get_coordinates_from_name(src_name2)
        assert cor1 == cor2
        assert cor1.to_string('hmsdms') == '00h00m41.527583s +39d18m04.14836s'
        src1 = sources.Source.source_from_name(src_name1)
        src2 = sources.Source.source_from_name(src_name2)
        assert src1.coord == src2.coord
        assert src1.coord == cor2
    except coord.name_resolve.NameResolveError:
        rprint("[yellow]WARNING: seems like no internet connection is available. "
               "Unable to resolve source names[/yellow]")


def test_scan_block():
    target_source = sources.Source('target', '10h20m10s 40d30m10s', source_type=sources.SourceType.TARGET)
    pcal1_source = sources.Source('pcal1', '10h21m00s 40d31m00s', source_type=sources.SourceType.PHASECAL)
    pcal2_source = sources.Source('pcal2', '10h19m00s 40d29m00s', source_type=sources.SourceType.PHASECAL)
    check_source = sources.Source('check', '10h20m00s 40d31m00s', source_type=sources.SourceType.CHECKSOURCE)
    target_scan = sources.Scan(source=target_source, duration=3.5*u.min)
    pcal1_scan = sources.Scan(source=pcal1_source, duration=1.5*u.min)
    pcal2_scan = sources.Scan(source=pcal2_source, duration=1.5*u.min)
    check_scan = sources.Scan(source=check_source, duration=3*u.min, every=3)
    with pytest.raises(ValueError):
        sources.ScanBlock([])
        sources.ScanBlock(['lfla', 2])

    block = sources.ScanBlock([pcal1_scan, target_scan, pcal2_scan, check_scan])
    assert len(block.scans) == 4
    assert block.has(sources.SourceType.TARGET)
    assert block.has(sources.SourceType.PHASECAL)
    assert block.has(sources.SourceType.CHECKSOURCE)
    assert len(block.sources(sources.SourceType.TARGET)) == 1
    assert len(block.sources(sources.SourceType.PHASECAL)) == 2
    assert len(block.sources(sources.SourceType.CHECKSOURCE)) == 1
    assert len(block.sources(sources.SourceType.FRINGEFINDER)) == 0

    scans = block.fill(1*u.h)
    assert len(scans) == 29
    assert reduce(operator.add, [s.duration for s in scans]) <= 1*u.h
    assert all([scan.source.type is sources.SourceType.CHECKSOURCE for scan in scans[8::9]])
