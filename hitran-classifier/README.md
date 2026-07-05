# HITRAN Classifier

Closes the gap between what JWST measures and what HITRAN can tell you about it.

JWST produces transmission and emission spectra of exoplanet atmospheres —
spectral fingerprints of every molecule present. HITRAN (the High Resolution
Transmission Molecular Absorption database) is the reference library that
tells you which molecule produced which feature. Matching one to the other
by hand, feature by feature, band by band, is slow, technical, and easy to
get wrong.

This tool automates that cross-referencing: give it a spectrum, it gives you
back ranked molecule candidates, flagged by biosignature significance.

It doesn't tell you there's life on a planet. Nothing can do that from a
spectrum alone. It makes the question *is this feature biological* faster
and more systematic to ask.

## What it does

1. **Loads** a spectrum (wavelength vs. transmission/flux/transit depth).
2. **Detects features** — absorption/emission peaks above the noise floor,
   using a rolling-baseline continuum estimate and prominence-based peak
   detection.
3. **Matches** each feature against a molecular reference table of band
   centers, scored by proximity.
4. **Ranks** candidate molecules and flags which ones are biosignature
   candidates — with a short note on *why* (biogenic vs. abiotic pathways,
   known false-positive risks, whether it needs a co-occurring gas to mean
   anything).

## Quick start

```bash
git clone https://github.com/<your-username>/hitran-classifier.git
cd hitran-classifier
pip install -r requirements.txt

# generates a synthetic demo spectrum with CH4/H2O/CO2 features baked in
python3 examples/generate_example_spectrum.py

# run the classifier on it
python3 -m hitran_classifier.cli examples/example_spectrum.csv

# with a plot
python3 -m hitran_classifier.cli examples/example_spectrum.csv --plot
```

## Using your own spectrum

Any CSV with a wavelength column (microns) and a signal column
(transit depth, transmission, or flux) works. Common column name
variants are auto-detected (`wavelength_um`, `wavelength`, `transit_depth`,
`flux`, etc). An optional error/uncertainty column improves feature
detection by setting the noise floor correctly.

```bash
python3 -m hitran_classifier.cli path/to/your_spectrum.csv
```

## Using a real HITRAN line list

The classifier ships with an approximate, curated reference table
(`data/biosignature_reference.json`) covering the molecules most relevant
to habitability and biosignature work — H2O, CO2, CH4, O2, O3, N2O,
dimethyl sulfide, methyl chloride, NH3, CO, H2S — with rough band centers
so it works offline, out of the box.

For rigorous line-by-line identification, plug in a real HITRAN export
instead (requires a free account at [hitran.org](https://hitran.org), or
use the [HITRAN Application Programming Interface](https://hitran.org/hapi/) — "hapi" — to pull line lists):

```bash
python3 -m hitran_classifier.cli examples/example_spectrum.csv \
    --hitran-file my_hitran_export.csv \
    --min-intensity 1e-24
```

The HITRAN export just needs molecule, and either wavenumber (cm⁻¹) or
wavelength (µm), columns — the loader collapses individual lines into
representative band centers automatically (see `hitran_classifier/hitran_loader.py`).

## Project layout

```
hitran_classifier/
  spectrum.py        # load + represent observed spectra
  features.py         # peak/feature detection above the noise floor
  matcher.py           # molecule matching, scoring, ranking
  hitran_loader.py    # optional: ingest real HITRAN line lists
  cli.py               # command-line interface
data/
  biosignature_reference.json   # built-in molecular reference table
examples/
  generate_example_spectrum.py  # synthetic demo spectrum generator
tests/
  test_matcher.py
```

## Why this exists

Two extraordinary tools,JWST and HITRAN,don't speak to each other
automatically. Matching an observed spectral feature to a molecular
identification means knowing which reference entries are relevant to your
wavelength range, your planet's temperature, your star's type.. work that's
manageable for a funded research team and a real barrier for anyone else
trying to engage seriously with biosignature data. This is infrastructure
to close that gap.

## Disclaimer

The built-in reference table uses approximate band centers for
demonstration and first-pass screening. It is **not** a substitute for
real line-by-line radiative transfer analysis. Any biosignature candidate
flagged here needs confirmation against the full HITRAN database (or
equivalent), atmospheric retrieval modeling, and for anything claiming
biological origin, a lot more scrutiny than one script can provide.

## License

MIT — see [LICENSE](LICENSE).
