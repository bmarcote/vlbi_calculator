"""Defines some of the possible values for the setup in a VLBI observation
"""
from typing import Union, Optional
from astropy import units as u

bands: dict[str, str] = {'92cm': '92 cm or 0.33 GHz', '49cm': '49 cm or 0.6 GHz',
                         '30cm': '30 cm or  1 GHz', '21cm': '21 cm or 1.4 GHz',
                         '18cm': '18 cm or 1.7 GHz', '13cm': '13 cm or 2.3 GHz',
                         '6cm': '6 cm or 5 GHz', '5cm': '5 cm or 6 GHz',
                         '3.6cm': '3.6 cm or 8.3 GHz', '2.5cm': '2.5 cm or 12 GHz',
                         '2cm': '2 cm or 15 GHz)',
                         '1.3cm': '1.3 cm or 23 GHz',  # '0.9cm': 'Ka band (0.9 cm or 33 GHz)',
                         '0.7cm': '0.7 cm or 43 GHz',
                         '0.348cm': '0.35 cm or 86 GHz',
                         '0.13cm': '0.13 cm or 230 GHz',
                         '0.09cm': '0.087 cm or 345 GHz'}

# from 4 Mbps to 32 Gbps
data_rates: dict[int, str] = {i: f"{i} Mbps" if i < 1e3 else f"{i/1000:.0f} Gbps"
                              for i in [2**j for j in range(2, 16)]}

# From 1 to 32 subbands (i.e. SPWs or IFs)
subbands: dict[int, str] = {2**i: f"{2**i} subbands" for i in range(6, -1, -1)}

# From 16 to 8192 spectral channels
channels: dict[int, str] = {2**i: f"{2**i} channels per subband" for i in range(14, 4, -1)}

# Full polarization meaning LL, RR, RL, LR; dual being RR and LL, and single just one of LL or RR.
polarizations: dict[int, str] = {4: '4 (Full polarization)', 2: '2 (dual polarization)',
                                 1: '1 (single polarization)'}

# Integration (averaging) times
inttimes: dict[Union[int, float], str] = {16: '16 s', 8: '8 s', 4: '4 s', 2: '2 s', 1: '1 s',
                                          0.5: '0.5 s', 0.25: '0.25 s', 0.125: '125 ms',
                                          0.06: '60 ms', 0.03: '30 ms', 0.015: '15 ms',
                                          0.001: '1 ms'}


def phaseref_cycle(band: str) -> Optional[u.Quantity]:
    """Returns the preferred phase referencing cycle for the given band for a standard VLBI observation.

    Input
        band :  str
            The observing band, e.g. '92cm'.
    Returns
        u.Quantity or None
            The preferred phase referencing cycle for the given band, e.g. 5*u.min at 6cm.
            Note that this cycle length means the time interval between two phase referencing scans,
            and thus it should be interpreted as the recommended length of a phasecal + target scans.
            It will return None if the observing frequency is > 60 GHz (wavelegnth < 0.5 cm), meaning
            no phase referencing is possible.
    """
    wavelength = float(band.replace('cm', ''))
    if wavelength > 15:
        return 7*u.min
    elif wavelength > 2:
        return 5*u.min
    elif wavelength > 0.5:
        return 2*u.min
    else:
        return None
