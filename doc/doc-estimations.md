


#### Data rate

VLBI observations are typically limited by the data rate: the amount of data that a station can transfer/record simultaneously. For simplification, this tool will always assume that all antennas record at the same speed. Whereas this is in general true, in multiple observations a few antennas can be limited to lower data rates. While this is fully supported in the SFXC correlator, the final effect is that those stations would show a fewer bandwidth, and thus the expected thermal rms noise may be slightly larger than expected.

The data rate for a given station is determined as
![equation2]({src:eq-datarate.png})
where _&#916;&#957;_ is the total observed bandwidth, _Np_ is the number of recorded polarizations (one or two), _Nb_ is the number of bits used to sample the data (VLBI observations typically record at 2-bit sampling), and the last 2 is related to the Nyquist sampling.



#### Thermal noise level

This tool assumes the nominal System Equivalent Flux Densities (SEFDs) for each antenna. The SEFD is define as the flux density of a radio source that doubles the system temperature and thus allows a direct conversion between the amplitudes recorded at each antenna and the real flux density of the emission recorded. The SEFD values for EVN antennas can be seen at [the EVN Status Table](http://old.evlbi.org/user_guide/EVNstatus.txt).

The thermal noise for a given observation can be estimated (following e.g. [Wrobel & Walker 1999](https://ui.adsabs.harvard.edu/abs/1999ASPC..180..171W/abstract)) by
![equation]({src:eq-noise.png})
where _&#414;_ is the system efficiency (typically 0.7 in VLBI, 2-bit, observations), _&#916;&#957;_ is the total bandwidth, _Np_ is the number of polarizations recorded at each antenna (one or two), _&#916;t_ is the total on-source integration time, and _SEFD_ is the SEFD for the _k_ station, where it is considered that there are _N_ antennas in total. Note that this tool takes into account that a given baseline (formed by the _i,j_ antennas) can observe for a limited amount of time, not necessarily the full observing time.

We remark that nominal SEFD values are considered. This implies that no gain-curve corrections are applied. This can have a significant effect when observing at low elevations, specially at higher frequencies. The SEFDs at 3 mm assume typical weather conditions and a mean elevation of 30 degrees. We also do not take into account expected losses due to external factors like radio frequency interferences (RFI; that can be significant at low frequencies). Therefore, the final rms noise level may result to be slightly larger than the thermal noise level.


