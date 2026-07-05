"""
HITRAN Classifier
=================

Connects JWST-style transmission/emission spectra to molecular
absorption/emission signatures, ranks candidate molecules by match
quality, and flags biosignature significance.

Built for the search for life elsewhere -- one gap closed at a time.
"""

__version__ = "0.1.0"

from .spectrum import load_spectrum, Spectrum
from .features import detect_features
from .matcher import classify_spectrum, MoleculeMatch
from .hitran_loader import load_hitran_linelist

__all__ = [
    "load_spectrum",
    "Spectrum",
    "detect_features",
    "classify_spectrum",
    "MoleculeMatch",
    "load_hitran_linelist",
]
