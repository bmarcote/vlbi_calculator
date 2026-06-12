"""vlbiplanobs package.

All public names are re-exported lazily (PEP 562) so that importing the package
(e.g. by the `planobs` console script) does not pull in the heavy dependencies
(astropy, numpy, catalogs) until they are actually used.
"""

__all__ = ["VLBIObs", "Time", "Stations", "Station",
           "Source", "Scan", "ScanBlock", "SourceType", "SourceNotVisible",
           "NETWORKS", "STATIONS",
           "RFCCatalog", "CalibratorSource", "get_fringe_finder_sources", "get_nearby_sources"]

# Maps each public name to (module, attribute) for lazy resolution.
_LAZY_ATTRS = {
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
        from importlib import import_module
        module_name, attr_name = _LAZY_ATTRS[name]
        value = getattr(import_module(module_name), attr_name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(__all__))
