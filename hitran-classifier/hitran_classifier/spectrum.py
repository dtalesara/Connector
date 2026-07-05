"""
Load and represent observed spectra (e.g. JWST transmission spectra).

Expected input: a CSV with at minimum a wavelength column (microns) and
a signal column (transit depth, transmission, or flux). Column names
are flexible -- common aliases are auto-detected.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd

WAVELENGTH_ALIASES = ["wavelength", "wavelength_um", "wave", "lambda", "um", "micron", "microns"]
SIGNAL_ALIASES = [
    "transit_depth", "depth", "transmission", "flux", "signal",
    "rp2_rs2", "(rp/rs)^2", "value", "y",
]
ERROR_ALIASES = ["error", "err", "sigma", "uncertainty", "yerr"]


@dataclass
class Spectrum:
    wavelength_um: np.ndarray
    signal: np.ndarray
    error: Optional[np.ndarray] = None
    label: str = "spectrum"

    def __post_init__(self):
        order = np.argsort(self.wavelength_um)
        self.wavelength_um = np.asarray(self.wavelength_um)[order]
        self.signal = np.asarray(self.signal)[order]
        if self.error is not None:
            self.error = np.asarray(self.error)[order]

    def __len__(self):
        return len(self.wavelength_um)


def _find_column(columns, aliases):
    lowered = {c.lower().strip(): c for c in columns}
    for alias in aliases:
        if alias in lowered:
            return lowered[alias]
    # fall back: substring match
    for alias in aliases:
        for lc, orig in lowered.items():
            if alias in lc:
                return orig
    return None


def load_spectrum(path: str, label: Optional[str] = None) -> Spectrum:
    """
    Load a spectrum from a CSV file.

    The loader tries to auto-detect wavelength / signal / error columns
    by name. If detection fails, it falls back to assuming the first
    column is wavelength and the second is signal.
    """
    df = pd.read_csv(path)
    wl_col = _find_column(df.columns, WAVELENGTH_ALIASES)
    sig_col = _find_column(df.columns, SIGNAL_ALIASES)
    err_col = _find_column(df.columns, ERROR_ALIASES)

    if wl_col is None or sig_col is None:
        if df.shape[1] < 2:
            raise ValueError(
                f"Could not identify wavelength/signal columns in {path}, "
                f"and file has fewer than 2 columns."
            )
        wl_col = wl_col or df.columns[0]
        sig_col = sig_col or df.columns[1]

    error = df[err_col].to_numpy(dtype=float) if err_col else None

    return Spectrum(
        wavelength_um=df[wl_col].to_numpy(dtype=float),
        signal=df[sig_col].to_numpy(dtype=float),
        error=error,
        label=label or path,
    )
