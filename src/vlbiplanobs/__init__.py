"""vlbiplanobs package.

All public names are re-exported lazily (PEP 562) so that importing the package
(e.g. by the `planobs` console script) does not pull in the heavy dependencies
(astropy, numpy, catalogs) until they are actually used.
"""
from importlib import import_module
__all__: list[str] = ["VLBIObs", "Time", "Stations", "Station",
           "Source", "Scan", "ScanBlock", "SourceType", "SourceNotVisible",
           "NETWORKS", "STATIONS",
           "RFCCatalog", "CalibratorSource", "get_fringe_finder_sources", "get_nearby_sources"]

# Maps each public name to (module, attribute) for lazy resolution.
_LAZY_ATTRS: dict[str, tuple[str, str]] = {
    "VLBIObs": ("vlbiplanobs.cli", "main"),
    "Time": ("astropy.time", "Time"),
    "Stations": ("vlbiplanobs.stations", "Stations"),
    "Station": ("vlbiplanobs.stations", "Station"),
    "Source": ("vlbiplanobs.sources", "Source"),
    "Scan": ("vlbiplanobs.sources", "Scan"),
    "ScanBlock": ("vlbiplanobs.sources", "ScanBlock"),
    "SourceType": ("vlbiplanobs.sources", "SourceType"),
    "SourceNotVisible": ("vlbiplanobs.sources", "SourceNotVisible"),
    "NETWORKS": ("vlbiplanobs.observation", "_NETWORKS"),
    "STATIONS": ("vlbiplanobs.observation", "_STATIONS"),
    "RFCCatalog": ("vlbiplanobs.calibrators", "RFCCatalog"),
    "CalibratorSource": ("vlbiplanobs.calibrators", "CalibratorSource"),
    "get_fringe_finder_sources": ("vlbiplanobs.calibrators", "get_fringe_finder_sources"),
    "get_nearby_sources": ("vlbiplanobs.calibrators", "get_nearby_sources"),
}


def __getattr__(name: str):
    """Lazily import and return the requested public attribute (PEP 562)."""
    if name in _LAZY_ATTRS:
        globals()[name] = getattr(import_module(name=_LAZY_ATTRS[name][0]), _LAZY_ATTRS[name][1])
        return globals()[name]

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(__all__))
