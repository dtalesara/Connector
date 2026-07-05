"""
Optional loader for real HITRAN line-list data.

The default biosignature_reference.json ships with approximate band
centers so the classifier works offline, out of the box. If you have
access to real HITRAN data (via an account at hitran.org, or exported
through the HITRAN Application Programming Interface, "hapi"), this
module lets you build a proper line-position reference table from it
and use it in place of the approximate defaults.

Supported input: a CSV/TSV export with at least these columns
(case-insensitive, common HITRAN/hapi export names are recognized):
    molecule / formula / molec_id
    wavenumber (cm^-1)  OR  wavelength_um
    intensity / line_intensity (optional, used for filtering weak lines)

If your export only has wavenumber, it is converted to microns via
wavelength_um = 1e4 / wavenumber_cm-1.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import pandas as pd

MOLECULE_COL_ALIASES = ["molecule", "formula", "molec_id", "species", "name"]
WAVENUMBER_COL_ALIASES = ["wavenumber", "nu", "wn", "cm-1", "cm^-1"]
WAVELENGTH_COL_ALIASES = ["wavelength_um", "wavelength", "lambda_um", "microns"]
INTENSITY_COL_ALIASES = ["intensity", "line_intensity", "sw", "s"]


def _find_col(columns, aliases):
    lowered = {c.lower().strip(): c for c in columns}
    for alias in aliases:
        if alias in lowered:
            return lowered[alias]
    for alias in aliases:
        for lc, orig in lowered.items():
            if alias in lc:
                return orig
    return None


def load_hitran_linelist(
    path: str,
    min_intensity: Optional[float] = None,
    band_tolerance_um: float = 0.02,
) -> Dict[str, List[float]]:
    """
    Load a real HITRAN line list export and collapse it into per-molecule
    band centers (microns), suitable for feeding into the matcher as a
    replacement reference table.

    Returns a dict: {molecule_name: [band_center_um, ...]}

    This intentionally simplifies a full line-by-line database into
    representative band centers, mirroring how the built-in reference
    table is structured, so it's a drop-in upgrade path rather than a
    parallel code path.
    """
    df = pd.read_csv(path, sep=None, engine="python")

    mol_col = _find_col(df.columns, MOLECULE_COL_ALIASES)
    wn_col = _find_col(df.columns, WAVENUMBER_COL_ALIASES)
    wl_col = _find_col(df.columns, WAVELENGTH_COL_ALIASES)
    inten_col = _find_col(df.columns, INTENSITY_COL_ALIASES)

    if mol_col is None:
        raise ValueError(
            "Could not find a molecule/formula column in the HITRAN export. "
            "Expected one of: " + ", ".join(MOLECULE_COL_ALIASES)
        )
    if wl_col is None and wn_col is None:
        raise ValueError(
            "Could not find a wavelength or wavenumber column in the HITRAN export."
        )

    if wl_col is None:
        df["_wavelength_um"] = 1e4 / df[wn_col].astype(float)
        wl_col = "_wavelength_um"

    if min_intensity is not None and inten_col is not None:
        df = df[df[inten_col].astype(float) >= min_intensity]

    bands: Dict[str, List[float]] = {}
    for molecule, group in df.groupby(mol_col):
        wavelengths = sorted(group[wl_col].astype(float).tolist())
        collapsed: List[float] = []
        for w in wavelengths:
            if not collapsed or (w - collapsed[-1]) > band_tolerance_um:
                collapsed.append(w)
        bands[str(molecule)] = collapsed

    return bands
