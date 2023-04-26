#!/usr/bin/env python3

import enum



# class bands():
#     """Observing bands available for the array.
#     """
#     @property
#
#     def __init__(self, fromfile='../data/frequency_setups.py'):

bands = {'92cm': '92 cm or 0.33 GHz', '49cm': '49 cm or 0.6 GHz',
         '30cm': '30 cm or  1 GHz', '21cm': '21 cm or 1.4 GHz',
         '18cm': '18 cm or 1.7 GHz', '13cm': '13 cm or 2.3 GHz',
         '6cm': '6 cm or 5 GHz', '5cm': '5 cm or 6 GHz',
         '3.6cm': '3.6 cm or 8.3 GHz', '2.5cm': '2.5 cm or 12 GHz',
         '2cm': '2 cm or 15 GHz)',
         '1.3cm': '1.3 cm or 23 GHz', #'0.9cm': 'Ka band (0.9 cm or 33 GHz)',
         '0.7cm': '0.7 cm or 43 GHz', '0.3cm': '0.3 cm or 100 GHz',
         '0.13cm': '0.13 cm or 230 GHz',
         '0.09cm': '0.087 cm or 345 GHz'}

# from 4 Mbps to 32 Gbps
data_rates = {i: f"{i} Mbps" if i<1e3 else f"{i/1000:.0f} Gbps" for i in [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 3000, 4096, 8192, 16384, 32768]}
subbands = {2**i: f"{2**i} subbands" for i in range(5, -1, -1)} # from 1 to 16
channels = {2**i: f"{2**i} channels per subband" for i in range(14, 4, -1)}# from 16 to 8192.
polarizations = {4: '4 (Full polarization)', 2: '2 (dual polarization)',
                 1: '1 (single polarization)'}
inttimes = {16: '16 s', 8: '8 s', 4: '4 s', 2: '2 s', 1: '1 s', 0.5: '0.5 s',
            0.25: '0.25 s', 0.125: '125 ms', 0.06: '60 ms', 0.03: '30 ms', 0.015: '15 ms',
            0.001: '1 ms'}



