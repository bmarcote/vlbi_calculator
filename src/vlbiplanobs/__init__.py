from .cli import main as VLBIObs
from astropy.time import Time
from .stations import Stations, Station
from .sources import Source, Scan, ScanBlock, SourceType, SourceNotVisible
from .observation import _NETWORKS as NETWORKS
from .observation import _STATIONS as STATIONS
