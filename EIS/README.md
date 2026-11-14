# EIS Fitting Module

Electrochemical Impedance Spectroscopy (EIS) equivalent circuit fitting with Constant Phase Elements (CPE) and Warburg diffusion elements.

## Project Objective

This module provides tools for fitting experimental EIS data to various equivalent circuit models, enabling extraction of electrochemical parameters such as solution resistance, charge transfer resistance, double-layer capacitance, and diffusion coefficients. The implementation includes multiple circuit topologies commonly encountered in corrosion science and electrochemistry.

## Technical Equations & Modeling

### Impedance Fundamentals

Complex impedance is defined as:
$$Z = Z' - jZ''$$

where $Z'$ is the real part (resistive component) and $Z''$ is the imaginary part (reactive component).

### Constant Phase Element (CPE)

The CPE impedance accounts for non-ideal capacitive behavior due to surface heterogeneity:

$$Z_{\text{CPE}} = \frac{1}{Q(j\omega)^n}$$

where:
- $Q$ is the CPE constant (s$^n$ / Ohm·cm²)
- $n$ is the dispersion coefficient (0 ≤ n ≤ 1)
  - $n = 1$: Ideal capacitor
  - $n = 0.5$: Warburg diffusion
  - $n < 1$: Non-ideal behavior

### Warburg Diffusion Elements

**Semi-infinite Warburg (infinite diffusion):**
$$Z_W = \frac{\sigma}{\sqrt{j\omega}}$$

where $\sigma$ is the Warburg coefficient related to diffusion coefficient.

**Finite-length Warburg (bounded diffusion):**
$$Z_W = R_W \frac{\tanh[(j\omega\tau_W)^{n_W}]}{(j\omega\tau_W)^{n_W}}$$

where:
- $R_W$ is the Warburg resistance
- $\tau_W$ is the diffusion time constant
- $n_W$ is typically 0.5 for Fickian diffusion

### Equivalent Circuit Models

#### Model L: Simplest Circuit
$$Z = R_s + \frac{1}{Q_1(j\omega)^{n_1} + \frac{1}{R_1}}$$

#### Model F: Series-Parallel
$$Z = R_s + \frac{1}{Q_1(j\omega)^{n_1} + \frac{1}{R_1}} + \frac{1}{Q_2(j\omega)^{n_2} + \frac{1}{R_2}}$$

#### Model J: Nested Parallel
$$Z = R_s + \frac{1}{Q_1(j\omega)^{n_1} + \frac{1}{R_1 + \frac{1}{Q_2(j\omega)^{n_2} + \frac{1}{R_2}}}}$$

#### Model B: Nested with Warburg
$$Z = R_s + \frac{1}{Q_1(j\omega)^{n_1} + \frac{1}{R_1 + \frac{1}{Q_2(j\omega)^{n_2} + \frac{1}{R_2 + Z_W}}}}$$

### Fitting Methodology

The fitting uses weighted least-squares optimization with impedance magnitude weighting:

$$\chi^2 = \sum_i \left[\frac{(Z'_{\text{exp},i} - Z'_{\text{model},i})^2 + (Z''_{\text{exp},i} - Z''_{\text{model},i})^2}{|Z_{\text{exp},i}|^2}\right]$$

Multi-start optimization with varying CPE exponents helps escape local minima.

## Usage

### Basic Example

```python
from eis_fitting import load_eis_data, fit_eis_model, z_model_series_parallel, plot_nyquist
import numpy as np

# Load data
freqs, Zr, Zi = load_eis_data('example_eis.csv')
omega = 2 * np.pi * freqs
Z_exp = Zr - 1j * Zi

# Define parameters
p0 = [10.0, 1e-6, 0.85, 1000.0, 1e-5, 0.75, 5000.0]
bounds = ([1.0, 1e-9, 0.3, 1.0, 1e-9, 0.3, 1.0],
          [100.0, 1e-2, 1.0, 1e6, 1e-2, 1.0, 1e7])

# Fit Model F (Series-Parallel)
popt, pct_err, chi2, chi2_red = fit_eis_model(
    z_model_series_parallel, p0, bounds, omega, Z_exp, lin_indices={2, 5}
)

# Plot results
Z_fit = z_model_series_parallel(omega, *popt)
plot_nyquist(Zr, Zi, Z_fit, freqs, title="Nyquist Plot - Model F Fit")
```

### Available Circuit Models

- `z_model_cpe_parallel`: Rs - (CPE1 // R1)
- `z_model_warburg_parallel`: Rs - CPE1 // (R1 - W)
- `z_model_series_parallel`: Rs - (CPE1 // R1) - (CPE2 // R2)
- `z_model_nested_parallel`: Rs - (CPE1 // (R1 + (CPE2 // R2)))
- `z_model_nested_diffusion`: Rs - CPE1 // (R1 - CPE2 // (R2 - W))
- `z_model_finite_warburg`: Rs - CPE1 // (R1 - CPE2 // (R2 - finite W))
- `z_model_triple_parallel`: Rs - CPE1 // (R1 - (R2 // CPE2) - (R3 // CPE3))
- `z_model_triple_series`: Rs + (CPE1 // R1) + (CPE2 // R2) + (CPE3 // R3)

## Output

![Nyquist Plot](../assets/eis_nyquist_example.png)

![Bode Plot](../assets/eis_bode_example.png)

## Data Format

Expected CSV columns:
- `Frequency (Hz)`: Frequency in Hz
- `Re(Z) normalized (Ohm.cm²)`: Real impedance
- `-Im(Z) normalized (Ohm.cm²)`: Negative imaginary impedance
- `Cycle` (optional): Cycle number for multi-cycle experiments

## Notes

- All impedances are area-normalized (Ohm·cm²)
- CPE exponents are bounded between 0 and 1
- Multi-start optimization improves convergence
- Error propagation via Jacobian covariance matrix

