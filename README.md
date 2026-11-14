# Materials Characterization & Analysis 

This repository showcases the scientific programming, data analysis pipelines, and visualization tools developed during my tenure as a Project Trainee in the Materials Science Division at Bhabha Atomic Research Centre (BARC).

It provides a modular suite of Python scripts for automated signal conditioning, non-linear curve fitting, and publication-quality visualization across electrochemical and surface-characterization techniques used in high-temperature oxidation and passive film research.

---

## Repository Architecture

| Module | Technique | Core Capabilities & Modeling Methods |
| :--- | :--- | :--- |
| **[`EIS/`](./EIS)** | Electrochemical Impedance Spectroscopy | Equivalent circuit fitting using Constant Phase Elements (CPE), Bode/Nyquist visualization. |
| **[`GDOES/`](./GDOES)** | Glow Discharge Optical Emission Spectroscopy | Sputter-rate depth calibration, compositional normalization, and elemental depth-profile visualization for oxide scale thickness and metal/oxide interface mapping. |
| **[`Raman/`](./Raman)** | Raman Spectroscopy | Asymmetric Least Squares (ALS) fluorescence baseline subtraction, Savitzky-Golay smoothing, and constrained multi-peak Pseudo-Voigt deconvolution with integrated phase area fractions. |
| **[`XRD/`](./XRD)** | X-Ray Diffraction | Automated multi-duration diffractogram alignment, dynamic vertical waterfall stacking, journal-formatted axes, and crystallographic peak indexing annotations. |

