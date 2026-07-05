"""
Command-line entry point.

Usage:
    python -m hitran_classifier.cli examples/example_spectrum.csv
    python -m hitran_classifier.cli examples/example_spectrum.csv --hitran-file my_hitran_export.csv
    python -m hitran_classifier.cli examples/example_spectrum.csv --plot
"""

from __future__ import annotations

import argparse
import sys

from .spectrum import load_spectrum
from .features import detect_features
from .matcher import classify_spectrum
from .hitran_loader import load_hitran_linelist


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="hitran-classifier",
        description="Match JWST-style spectra against molecular biosignature bands.",
    )
    p.add_argument("spectrum_csv", help="Path to a CSV file with wavelength + signal columns.")
    p.add_argument(
        "--hitran-file", default=None,
        help="Optional real HITRAN line-list export to use instead of/alongside the built-in reference table.",
    )
    p.add_argument(
        "--min-intensity", type=float, default=None,
        help="Minimum line intensity to keep when loading a HITRAN file (filters weak lines).",
    )
    p.add_argument(
        "--prominence", type=float, default=None,
        help="Minimum feature prominence for peak detection. Auto-set if omitted.",
    )
    p.add_argument("--plot", action="store_true", help="Show a matplotlib plot of the spectrum with matches.")
    p.add_argument("--top", type=int, default=10, help="Number of top matches to display.")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    spectrum = load_spectrum(args.spectrum_csv)
    features = detect_features(spectrum, prominence=args.prominence)

    hitran_bands = None
    if args.hitran_file:
        hitran_bands = load_hitran_linelist(args.hitran_file, min_intensity=args.min_intensity)

    result = classify_spectrum(spectrum, features=features, hitran_bands=hitran_bands)

    print(f"\nSpectrum: {result.spectrum_label}")
    print(f"Detected features: {len(result.features)}")
    for f in result.features:
        print(f"  {f.kind:10s}  {f.wavelength_um:.3f} um   strength={f.strength:.4f}")

    print(f"\nTop molecule matches (showing up to {args.top}):\n")
    header = f"{'Molecule':<22}{'Formula':<16}{'Score':>7}  {'Biosig':<7} {'Matched (um)'}"
    print(header)
    print("-" * len(header))
    for m in result.ranked()[: args.top]:
        biosig_flag = "YES" if m.biosignature else "no"
        wl_str = ", ".join(f"{w:.2f}" for w in m.matched_features_um)
        print(f"{m.molecule:<22}{m.formula:<16}{m.score:>7.2f}  {biosig_flag:<7} {wl_str}")

    biosigs = result.biosignature_candidates()
    if biosigs:
        print("\nBiosignature candidates worth a closer look:\n")
        for m in biosigs:
            print(f"- {m.molecule} ({m.formula}), score={m.score:.2f}")
            print(f"    {m.significance}")
    else:
        print("\nNo biosignature-flagged molecules matched above threshold.")

    if args.plot:
        _plot(spectrum, result)

    return 0


def _plot(spectrum, result):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(spectrum.wavelength_um, spectrum.signal, color="black", lw=1, label=spectrum.label)
    for feat in result.features:
        color = "crimson" if feat.kind == "absorption" else "steelblue"
        ax.axvline(feat.wavelength_um, color=color, alpha=0.3, lw=1)
    for m in result.ranked()[:5]:
        for wl in m.matched_features_um:
            ax.annotate(
                m.formula, xy=(wl, spectrum.signal.min()),
                xytext=(wl, spectrum.signal.min()),
                rotation=90, fontsize=8, ha="center",
            )
    ax.set_xlabel("Wavelength (microns)")
    ax.set_ylabel("Signal")
    ax.set_title(f"HITRAN Classifier — {spectrum.label}")
    ax.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    sys.exit(main())
