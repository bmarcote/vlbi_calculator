#!/usr/bin/env python3

import enum



# class bands():
#     """Observing bands available for the array.
#     """
#     @property
#
#     def __init__(self, fromfile='../data/frequency_setups.py'):

bands = {'92cm': 'P - 92 cm / 0.33 GHz', '49cm': 'P - 49 cm / 0.6 GHz',
         '30cm': 'UHF - 30 cm /  1 GHz', '21cm': 'L - 21 cm / 1.4 GHz',
         '18cm': 'L - 18 cm / 1.7 GHz', '13cm': 'S - 13 cm / 2.3 GHz',
         '6cm': 'C - 6cm / 5 GHz', '5cm': 'C - 5 cm / 6 GHz',
         '3.6cm': 'X - 3.6 cm / 8.3 GHz', '2cm': 'U - 2 cm / 15 GHz',
         '1.3cm': 'K - 1.3 cm / 23 GHz', '0.9cm': 'Ka - 0.9 cm / 33 GHz',
         '0.7cm': 'Q - 0.7 cm / 43 GHz', '0.3cm': 'W - 0.3 cm / 100 GHz',
         '0.1cm': '0.1 cm / 300 GHz'}

data_rates = (2**i for i in range(13, 2, -1)) # from 4 to 4096
subbands = (2**i for i in range(5, 0, -1)) # from 1 to 16
channels = (2**i for i in range(14, 4, -1)) # from 16 to 8192.
polarizations = {4: '4 (Full polarization)', 2: '2 (dual polarization)',
                 1: '1 (single polarization)'}
inttimes = {16: '16 s', 8: '8 s', 4: '4 s', 2: '2 s', 1: '1 s', 0.5: '0.5 s',
            0.25: '0.25 s', 0.125: '125 ms', 0.06: '60 ms', 0.03: '30 ms', 0.015: '15 ms',
            0.001: '1 ms'}



