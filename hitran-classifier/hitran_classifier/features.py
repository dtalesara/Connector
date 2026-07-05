"""
Detect candidate spectral features (absorption or emission peaks) in an
observed spectrum, ready to be matched against molecular reference bands.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List

import numpy as np
from scipy.signal import find_peaks

from .spectrum import Spectrum


@dataclass
class Feature:
    wavelength_um: float
    strength: float          # deviation from local baseline
    width_um: float          # rough feature width
    kind: str                # "absorption" or "emission"


def _rough_baseline(signal: np.ndarray, window: int) -> np.ndarray:
    """Simple rolling-median baseline for continuum estimation."""
    window = max(3, window | 1)  # ensure odd, >=3
    pad = window // 2
    padded = np.pad(signal, pad, mode="edge")
    baseline = np.array([
        np.median(padded[i:i + window]) for i in range(len(signal))
    ])
    return baseline


def detect_features(
    spectrum: Spectrum,
    prominence: float = None,
    baseline_window: int = 25,
    min_distance_points: int = 8,
) -> List[Feature]:
    """
    Detect absorption and emission features relative to a local
    continuum estimate.

    prominence: minimum prominence passed to scipy.find_peaks. If None,
    it's set automatically to 2.5 * the error array (if provided) or
    2.5 * std of the continuum-subtracted signal -- tuned to reject
    noise-level wiggles rather than flag every ripple as a "feature".
    min_distance_points: minimum spacing between accepted peaks, in
    samples, to avoid splitting one physical feature into many.
    """
    if len(spectrum) < 5:
        raise ValueError("Spectrum too short to detect features (need >= 5 points).")

    baseline = _rough_baseline(spectrum.signal, baseline_window)
    residual = spectrum.signal - baseline

    if prominence is None:
        if spectrum.error is not None and np.any(spectrum.error > 0):
            noise_level = float(np.median(spectrum.error))
        else:
            noise_level = float(np.std(residual))
        prominence = 2.5 * noise_level if noise_level > 0 else 1e-6

    features: List[Feature] = []

    # Absorption features: residual dips (negative for transmission-style
    # signals where absorption reduces flux; for transit depth spectra
    # absorption instead shows as a rise, both directions are checked)
    dip_idx, dip_props = find_peaks(-residual, prominence=prominence, distance=min_distance_points)
    for i, idx in enumerate(dip_idx):
        width = _estimate_width(spectrum.wavelength_um, idx, dip_props, i)
        features.append(Feature(
            wavelength_um=float(spectrum.wavelength_um[idx]),
            strength=float(abs(residual[idx])),
            width_um=width,
            kind="absorption",
        ))

    peak_idx, peak_props = find_peaks(residual, prominence=prominence, distance=min_distance_points)
    for i, idx in enumerate(peak_idx):
        width = _estimate_width(spectrum.wavelength_um, idx, peak_props, i)
        features.append(Feature(
            wavelength_um=float(spectrum.wavelength_um[idx]),
            strength=float(abs(residual[idx])),
            width_um=width,
            kind="emission",
        ))

    features.sort(key=lambda f: f.wavelength_um)
    return features


def _estimate_width(wavelength, idx, props, i, default=0.1):
    """Rough width estimate in microns based on array spacing near the peak."""
    n = len(wavelength)
    lo = max(0, idx - 3)
    hi = min(n - 1, idx + 3)
    if hi > lo:
        return float(wavelength[hi] - wavelength[lo]) / 2.0
    return default
