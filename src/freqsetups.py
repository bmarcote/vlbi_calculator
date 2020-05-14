#!/usr/bin/env python3

import enum



# class bands():
#     """Observing bands available for the array.
#     """
#     @property
#
#     def __init__(self, fromfile='../data/frequency_setups.py'):

bands = {'92cm': 'P - 92cm', '49cm': 'P - 49cm', '30cm': 'UHF - 30cm', '21cm': 'L - 21cm',
         '18cm': 'L - 18cm', '13cm': 'S - 13cm', '6cm': 'C - 6cm', '5cm': 'C - 5cm',
         '3.6cm': 'X - 3.6cm', '2cm': 'U - 2cm', '1.3cm': 'K - 1.3cm', '0.9cm': 'Ka - 0.9cm',
         '0.7cm': 'Q - 0.7cm', '0.3cm': 'W - 0.3cm', '0.1cm': '0.1cm'}

data_rates = (2**i for i in range(13, 2, -1)) # from 4 to 4096
subbands = (2**i for i in range(5, 0, -1)) # from 1 to 16
channels = (2**i for i in range(14, 4, -1)) # from 16 to 8192.
polarizations = {4: '4 (Full polarization)', 2: '2 (dual polarization)',
                 1: '1 (single polarization)'}
inttimes = {16: '16 s', 8: '8 s', 4: '4 s', 2: '2 s', 1: '1 s', 0.5: '0.5 s',
            0.25: '0.25 s', 0.125: '125 ms', 0.06: '60 ms', 0.03: '30 ms', 0.015: '15 ms',
            0.001: '1 ms'}



