"""
Core classifier: matches detected spectral features against a molecular
reference table (default: built-in biosignature reference; optionally: a
real HITRAN-derived table), scores candidates, and ranks them with a
biosignature significance flag.

This is the piece that closes the gap described in the show notes --
JWST produces spectra, HITRAN (or an approximation of it) tells you what
molecules produce which features, and this module does the cross-referencing
so you don't have to do it by hand for every feature.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .features import Feature
from .spectrum import Spectrum
from .features import detect_features

DEFAULT_REFERENCE_PATH = os.path.join(
    os.path.dirname(__file__), "..", "data", "biosignature_reference.json"
)


@dataclass
class MoleculeMatch:
    molecule: str
    formula: str
    matched_features_um: List[float]
    score: float                      # 0-1, higher = better match
    biosignature: bool
    abiotic_source: bool
    significance: str
    reference_note: str


@dataclass
class ClassificationResult:
    spectrum_label: str
    features: List[Feature]
    matches: List[MoleculeMatch] = field(default_factory=list)

    def ranked(self) -> List[MoleculeMatch]:
        return sorted(self.matches, key=lambda m: m.score, reverse=True)

    def biosignature_candidates(self) -> List[MoleculeMatch]:
        return [m for m in self.ranked() if m.biosignature and m.score > 0]


def _load_default_reference(path: str = DEFAULT_REFERENCE_PATH) -> list:
    with open(path, "r") as f:
        data = json.load(f)
    return data["molecules"]


def _score_match(feature_wl: float, band_centers: List[float], band_width: float) -> float:
    """
    Score how well a single feature matches a molecule's set of band
    centers. Uses a Gaussian-ish proximity score against the nearest band.
    """
    if not band_centers:
        return 0.0
    nearest = min(band_centers, key=lambda b: abs(b - feature_wl))
    distance = abs(nearest - feature_wl)
    if distance > 3 * band_width:
        return 0.0
    # Score decays smoothly from 1.0 (exact match) to 0.0 (3 widths away)
    return max(0.0, 1.0 - (distance / (3 * band_width)))


def classify_spectrum(
    spectrum: Spectrum,
    features: Optional[List[Feature]] = None,
    reference: Optional[List[dict]] = None,
    hitran_bands: Optional[Dict[str, List[float]]] = None,
    match_tolerance_um: Optional[float] = None,
) -> ClassificationResult:
    """
    Classify a spectrum against the molecular reference table.

    Args:
        spectrum: the observed Spectrum.
        features: pre-detected features; if None, detect_features() is run.
        reference: override for the built-in biosignature reference table
            (list of molecule dicts, same schema as biosignature_reference.json).
        hitran_bands: optional dict of {molecule: [band_centers_um]} from a
            real HITRAN export (see hitran_loader.load_hitran_linelist).
            When provided, these bands are matched using match_tolerance_um
            as the band width, and results are merged with the reference
            table's biosignature metadata where names align.
        match_tolerance_um: band half-width to use for hitran_bands matches.
            Defaults to 0.02 microns if not given.

    Returns:
        ClassificationResult with detected features and ranked molecule matches.
    """
    if features is None:
        features = detect_features(spectrum)

    if reference is None:
        reference = _load_default_reference()

    matches: List[MoleculeMatch] = []

    # 1. Match against the curated biosignature reference table
    for mol in reference:
        band_centers = mol["bands_um"]
        band_width = mol.get("band_width_um", 0.15)
        best_score = 0.0
        matched_wavelengths = []
        for feat in features:
            s = _score_match(feat.wavelength_um, band_centers, band_width)
            if s > 0:
                matched_wavelengths.append(feat.wavelength_um)
                best_score = max(best_score, s)
        if matched_wavelengths:
            matches.append(MoleculeMatch(
                molecule=mol["name"],
                formula=mol["formula"],
                matched_features_um=sorted(set(round(w, 3) for w in matched_wavelengths)),
                score=round(best_score, 3),
                biosignature=mol["biosignature"],
                abiotic_source=mol["abiotic_source"],
                significance=mol["significance"],
                reference_note=mol["reference"],
            ))

    # 2. Optionally also match against a real HITRAN-derived band table
    if hitran_bands:
        width = match_tolerance_um or 0.02
        ref_lookup = {m["name"].lower(): m for m in reference}
        ref_lookup.update({m["formula"].lower(): m for m in reference})
        for molecule, band_centers in hitran_bands.items():
            best_score = 0.0
            matched_wavelengths = []
            for feat in features:
                s = _score_match(feat.wavelength_um, band_centers, width)
                if s > 0:
                    matched_wavelengths.append(feat.wavelength_um)
                    best_score = max(best_score, s)
            if matched_wavelengths:
                meta = ref_lookup.get(molecule.lower())
                matches.append(MoleculeMatch(
                    molecule=molecule,
                    formula=meta["formula"] if meta else molecule,
                    matched_features_um=sorted(set(round(w, 3) for w in matched_wavelengths)),
                    score=round(best_score, 3),
                    biosignature=meta["biosignature"] if meta else False,
                    abiotic_source=meta["abiotic_source"] if meta else True,
                    significance=meta["significance"] if meta else
                        "Matched via real HITRAN line list; no curated biosignature metadata available for this species.",
                    reference_note=meta["reference"] if meta else "Source: user-supplied HITRAN export.",
                ))

    return ClassificationResult(
        spectrum_label=spectrum.label,
        features=features,
        matches=matches,
    )
