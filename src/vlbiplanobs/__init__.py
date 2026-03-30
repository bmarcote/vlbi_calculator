__version__ = "4.8.1"

from .cli import main as VLBIObs
from astropy.time import Time
from .stations import Stations, Station
from .sources import Source, Scan, ScanBlock, SourceType, SourceNotVisible
from .observation import _NETWORKS as NETWORKS
from .observation import _STATIONS as STATIONS
from .calibrators import (
    RFCCatalog,
    CalibratorSource,
    get_fringe_finder_sources,
    get_nearby_sources
)

__all__ = [
    "__version__",
    "VLBIObs", "Time", "Stations", "Station",
    "Source", "Scan", "ScanBlock", "SourceType", "SourceNotVisible",
    "NETWORKS", "STATIONS",
    "RFCCatalog", "CalibratorSource", "get_fringe_finder_sources", "get_nearby_sources",
]
