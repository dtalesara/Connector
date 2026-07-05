import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hitran_classifier.spectrum import Spectrum
from hitran_classifier.features import detect_features, Feature
from hitran_classifier.matcher import classify_spectrum


def make_synthetic_spectrum_with_ch4_dip():
    # Wide wavelength range with a narrow feature, so the rolling-median
    # baseline (window=25 samples) stays wider than the feature itself
    # and can actually estimate a continuum instead of tracking the dip.
    wavelength_um = np.linspace(2.5, 4.5, 800)
    baseline = np.full_like(wavelength_um, 0.01)
    dip = 0.002 * np.exp(-0.5 * ((wavelength_um - 3.3) / 0.015) ** 2)
    signal = baseline + dip
    error = np.full_like(wavelength_um, 0.00005)
    return Spectrum(wavelength_um=wavelength_um, signal=signal, error=error, label="test")


def test_detect_features_finds_the_injected_dip():
    spectrum = make_synthetic_spectrum_with_ch4_dip()
    features = detect_features(spectrum)
    assert len(features) >= 1
    closest = min(features, key=lambda f: abs(f.wavelength_um - 3.3))
    assert abs(closest.wavelength_um - 3.3) < 0.05


def test_classify_spectrum_flags_methane_as_biosignature_candidate():
    spectrum = make_synthetic_spectrum_with_ch4_dip()
    result = classify_spectrum(spectrum)
    names = [m.molecule for m in result.matches]
    assert "Methane" in names
    methane_match = next(m for m in result.matches if m.molecule == "Methane")
    assert methane_match.biosignature is True
    assert methane_match.score > 0.5


def test_no_features_means_no_matches():
    wavelength_um = np.linspace(3.0, 3.6, 200)
    flat_signal = np.full_like(wavelength_um, 0.01)
    spectrum = Spectrum(wavelength_um=wavelength_um, signal=flat_signal, label="flat")
    result = classify_spectrum(spectrum)
    assert len(result.features) == 0
    assert len(result.matches) == 0


def test_biosignature_candidates_only_returns_flagged_molecules():
    spectrum = make_synthetic_spectrum_with_ch4_dip()
    result = classify_spectrum(spectrum)
    for m in result.biosignature_candidates():
        assert m.biosignature is True


if __name__ == "__main__":
    test_detect_features_finds_the_injected_dip()
    test_classify_spectrum_flags_methane_as_biosignature_candidate()
    test_no_features_means_no_matches()
    test_biosignature_candidates_only_returns_flagged_molecules()
    print("All tests passed.")
