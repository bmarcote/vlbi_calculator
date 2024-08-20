# import os
# import sys
from dataclasses import dataclass
from typing import Self #, Optional, Union, Sequence
# import numpy as np
from astropy import coordinates as coord





@dataclass
class Source:
    name: str
    ivsname: str
    coordinates: coord.SkyCoord
    fluxes: dict[str, list[float]]
    category: str
    url: str




class CalibratorCatalog:
    """Full RFC catalog with all sources.
    """
    @property
    def sources(self) -> dict:
        """Returns all sources in the catalog.
        """
        raise NotImplementedError


    def unresolved(self) -> Self:
        """Returns all sources that remain unresolved at long baselines.
        """
        raise NotImplementedError


    def brighter_than(self, flux: float, band: str) -> Self:
        """Returns all sources that are brighter than the given amount on the
        longest baselines.

        Input
        - flux : float
            Flux density on the longest baselines in Jansky units.
        - band : str
            Band at which the flux is compared. It needs to be 'S', 'C', 'X', 'U', 'K'
        """
        raise NotImplementedError


    def get(self, source_name: str) -> Source:
        """Returns the source associated to the given name.
        """
        raise NotImplementedError


    def import_rfc(self):
        """Imports the RFC catalog of sources
        """
        raise NotImplementedError


    def import_vsop(self):
        """Imports the VSOP catalog of sources
        """
        raise NotImplementedError


    def import_icrf(self):
        """Imports the ICRF catalog of sources
        """
        raise NotImplementedError




