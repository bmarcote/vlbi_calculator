import timeit
import numpy as np
from vlbiplanobs import observation as obs
# from vlbiplanobs import stations
from astropy import units as u
from astropy.time import Time

def test_speed():
    o = obs.Observation()
    o.polarizations = 'dual'
    assert o.polarizations.value == 2
    o.polarizations = 4
    assert o.polarizations == obs.Polarization.FULL
    target = obs.Source('Target', '10h58m29.6s +81d33m58.8s')

    o = obs.Observation(target=target)
    o.times = Time('2020-06-15 20:00', scale='utc') + np.arange(0, 720, 10)*u.min
    o.band = '18'
    o.datarate = 1024
    o.subbands = 8
    o.channels = 32
    o.polarizations = 2
    o.inttime = 2

    # all_stations = stations.Network.get_stations_from_configfile()
    evn6 = ['Ef', 'Jb2', 'On', 'Hh', 'T6', 'Wb', 'Sv', 'Zc', 'Pa', 'Mp', 'Ho', 'Nl', 'Pt', 'Sc', 'Kp', 'Hn']
    o.stations_from_codenames(evn6)


    print(f"Elevations with process: {timeit.timeit(o._elevations_process, number=5)}")
    print(f"Elevations with threads: {timeit.timeit(o._elevations_threads, number=5)}")
    print(f"Elevations with nothing: {timeit.timeit(o._elevations_single, number=5)}")
    print(f"Elevations with asyncio: {timeit.timeit(o._elevations_asyncio, number=5)}")
    # Elevations with process: 15.168135124957189
    # Elevations with threads: 0.21708729199599475
    # Elevations with nothing: 0.36412112502148375
    # Elevations with asyncio: 0.36775329196825624
    # Altaz with process: 14.9550038339803
    # Altaz with threads: 0.20259183296002448
    # Altaz with nothing: 0.3541361659881659
    _ = o.elevations()

    print(f"Altaz with process: {timeit.timeit(o._altaz_process, number=5)}")
    print(f"Altaz with threads: {timeit.timeit(o._altaz_threads, number=5)}")
    print(f"Altaz with nothing: {timeit.timeit(o._altaz_single, number=5)}")
    _ = o.is_visible()
    # beam = o.synthesized_beam()
    # _ = o.thermal_noise()
    # _ = o.get_uv_array()
    # _ = o.get_dirtymap(robust='natural')
    # _ = o.get_dirtymap(robust='uniform')
    # fig, ax = plt.subplots()



if __name__ == '__main__':
    test_speed()
