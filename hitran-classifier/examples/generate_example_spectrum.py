"""
Generates a synthetic JWST-style transmission spectrum with absorption
dips placed near real CH4, H2O, and CO2 bands, plus noise -- so the
classifier has something realistic to chew on out of the box.

This is SYNTHETIC DATA for demonstration only, not a real observation.
"""

import numpy as np
import pandas as pd

rng = np.random.default_rng(42)

wavelength_um = np.linspace(1.0, 12.0, 600)
baseline = 0.0142 + 0.0002 * wavelength_um  # gentle continuum slope
signal = baseline.copy()

bands = [
    (1.4, 0.05, 0.0006),   # H2O
    (1.7, 0.05, 0.0004),   # CH4
    (2.3, 0.06, 0.0005),   # CH4
    (2.7, 0.07, 0.0009),   # H2O / CO2 blend
    (3.3, 0.06, 0.0007),   # CH4
    (4.3, 0.08, 0.0011),   # CO2
    (6.3, 0.10, 0.0005),   # H2O
    (7.7, 0.10, 0.0006),   # CH4
]

for center, width, depth in bands:
    signal += depth * np.exp(-0.5 * ((wavelength_um - center) / width) ** 2)

noise = rng.normal(0, 0.00012, size=wavelength_um.shape)
signal_noisy = signal + noise
error = np.full_like(wavelength_um, 0.00012)

df = pd.DataFrame({
    "wavelength_um": wavelength_um,
    "transit_depth": signal_noisy,
    "error": error,
})

df.to_csv("examples/example_spectrum.csv", index=False)
print("Wrote examples/example_spectrum.csv with", len(df), "rows")
