#!/usr/bin/env python3

# class bands():
#     """Observing bands available for the array.
#     """
#     @property
#
#     def __init__(self, fromfile='../data/frequency_setups.py'):

bands = {'92': '92 cm or 0.33 GHz', '49': '49 cm or 0.6 GHz',
         '30': '30 cm or  1 GHz', '21': '21 cm or 1.4 GHz',
         '18': '18 cm or 1.7 GHz', '13': '13 cm or 2.3 GHz',
         '6': '6 cm or 5 GHz', '5': '5 cm or 6 GHz',
         '3.6': '3.6 cm or 8.3 GHz', '2.5': '2.5 cm or 12 GHz',
         '2': '2 cm or 15 GHz)',
         '1.3': '1.3 cm or 23 GHz', #'0.9cm': 'Ka band (0.9 cm or 33 GHz)',
         '0.7': '0.7 cm or 43 GHz',
         '0.348': '0.35 cm or 86 GHz',
         '0.13': '0.13 cm or 230 GHz',
         '0.09': '0.087 cm or 345 GHz'}

# from 4 Mbps to 32 Gbps
data_rates = {i: f"{i} Mbps" if i<1e3 else f"{i/1000:.0f} Gbps" for i in [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 3000, 4096, 8192, 16384, 32768]}
subbands = {2**i: f"{2**i} subbands" for i in range(5, -1, -1)} # from 1 to 16
channels = {2**i: f"{2**i} channels per subband" for i in range(14, 4, -1)}# from 16 to 8192.
polarizations = {4: '4 (Full polarization)', 2: '2 (dual polarization)',
                 1: '1 (single polarization)'}
inttimes = {16: '16 s', 8: '8 s', 4: '4 s', 2: '2 s', 1: '1 s', 0.5: '0.5 s',
            0.25: '0.25 s', 0.125: '125 ms', 0.06: '60 ms', 0.03: '30 ms', 0.015: '15 ms',
            0.001: '1 ms'}



