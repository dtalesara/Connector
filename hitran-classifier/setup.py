from setuptools import setup, find_packages

setup(
    name="hitran-classifier",
    version="0.1.0",
    description="Match JWST-style exoplanet spectra against molecular biosignature bands.",
    packages=find_packages(exclude=["tests"]),
    install_requires=[
        "numpy>=1.24",
        "scipy>=1.10",
        "pandas>=2.0",
        "matplotlib>=3.7",
    ],
    python_requires=">=3.9",
    include_package_data=True,
    package_data={"": ["data/*.json"]},
    entry_points={
        "console_scripts": [
            "hitran-classifier=hitran_classifier.cli:main",
        ],
    },
)
