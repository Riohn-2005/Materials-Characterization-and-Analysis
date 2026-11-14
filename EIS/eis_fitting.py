
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from typing import Tuple, Dict, Callable, Optional

# Start by definig the circuit to be used
# Here I am using a simple circuit as an example.

# Model E
def Z_model_E(omega, Rs, Q1, n1, R1, sigma):
      """
    Simplest model: Rs - (CPE1 // R1 - W)
    Parameters:
        omega: Angular frequency (rad/s)
        Rs: Solution resistance (Ohm·cm²)
        Q1: CPE constant (s^n / Ohm·cm²)
        n1: CPE dispersion coefficient (0 ≤ n1 ≤ 1)
        R1: Polarization resistance (Ohm·cm²)
        sigma: Warburg
    """
    j = 1j
    Z_W = sigma / np.sqrt(j * omega)
    Z_branch = R1 + Z_W
    Z_outer = 1.0 / (Q1 * (j * omega)**n1 + 1.0 / Z_branch)
    return Rs + Z_outer
    

# This is an example of a complex circuit
# Model H: Nested Branch Circuit (Rs + (CPE1 // (R1 + CPE2 // (R2 + CPE3 // R3))))
def Z_model_H(omega, Rs, Q1, n1, R1, Q2, n2, R2, Q3, n3, R3):
    j = 1j
    Z_cpe3 = 1.0 / (Q3 * (j * omega)**n3)
    Z_inner2 = 1.0 / (1.0 / Z_cpe3 + 1.0 / R3)
    Z_cpe2 = 1.0 / (Q2 * (j * omega)**n2)
    Z_inner1 = 1.0 / (1.0 / Z_cpe2 + 1.0 / (R2 + Z_inner2))
    Z_cpe1 = 1.0 / (Q1 * (j * omega)**n1)
    Z_outer = 1.0 / (1.0 / Z_cpe1 + 1.0 / (R1 + Z_inner1))
    return Rs + Z_outer

# FITTING FUNCTIONS
# Then I create a function that takes the model to be fitted, intial guess, bounds (range) of my values and the experimental data
def fit_model(model_func, p0, bounds, omega, Z_exp):

    def residuals(log_params):
        Rs = np.exp(log_params[0])
        Q1 = np.exp(log_params[1])
        n1 = log_params[2]          # n1 is fitted linearly
        R1 = np.exp(log_params[3])
        sigma = np.exp(log_params[4])

        Zm = Z_model_E(omega, Rs, Q1, n1, R1, sigma)
        w = np.abs(Z_exp)
        return np.concatenate([(Z_exp.real - Zm.real)/w, (Z_exp.imag - Zm.imag)/w])

    # Convert initial guess to fitting space (log for magnitudes, linear for n)
    log_p0 = [np.log(p0[0]), np.log(p0[1]), p0[2], np.log(p0[3]), np.log(p0[4])]

    # Convert bounds to fitting space
    log_bounds = (
        [np.log(bounds[0][0]), np.log(bounds[0][1]), bounds[0][2], np.log(bounds[0][3]), np.log(bounds[0][4])],
        [np.log(bounds[1][0]), np.log(bounds[1][1]), bounds[1][2], np.log(bounds[1][3]), np.log(bounds[1][4])]
    )

    best_res = None
    best_chi2 = np.inf

    # Try multiple start points for the CPE exponent n1 to avoid local minima
    starts = [log_p0]
    for n1_v in [0.75, 0.85]:
        starts.append([log_p0[0], log_p0[1], n1_v, log_p0[3], log_p0[4]])

    for x0 in starts:
        x0 = np.clip(x0, log_bounds[0], log_bounds[1])
        try:
            res = least_squares(residuals, x0, bounds=log_bounds, method='trf', ftol=1e-15, xtol=1e-15, gtol=1e-15)
            chi2_val = np.sum(res.fun**2)
            if chi2_val < best_chi2:
                best_chi2 = chi2_val
                best_res = res
        except Exception as e:
            continue

    # Exponentiate the optimized log parameters back to real values
    log_opt = best_res.x
    popt = np.zeros(5)
    popt[0] = np.exp(log_opt[0])  # Rs
    popt[1] = np.exp(log_opt[1])  # Q1
    popt[2] = log_opt[2]          # n1
    popt[3] = np.exp(log_opt[3])  # R1
    popt[4] = np.exp(log_opt[4])  # sigma

    # Calculate statistics
    chi2 = np.sum(best_res.fun**2)
    n_points = len(omega)
    chi2_red = chi2 / (2 * n_points - 5)  # 5 degrees of freedom lost

    # Standard errors extraction via Jacobian
    J = best_res.jac
    dof = len(best_res.fun) - 5
    s2 = chi2 / max(dof, 1)
    cov_log = np.linalg.pinv(J.T @ J) * s2
    std_log = np.sqrt(np.maximum(0, np.diag(cov_log)))

    # Error propagation from log-space to linear-space
    perr = np.zeros(5)
    perr[0] = popt[0] * std_log[0]
    perr[1] = popt[1] * std_log[1]
    perr[2] = std_log[2]          # n1 standard error remains un-scaled
    perr[3] = popt[3] * std_log[3]
    perr[4] = popt[4] * std_log[4]

    pct_err = (perr / popt) * 100.0

    return popt, pct_err, chi2, chi2_red


# Base Estimates ---
Rs_est = Zr[0]
Rtotal_est = Zr[-1]
Rp_est = Rtotal_est - Rs_est

# Index Guide: [Rs, Q1, n1, R1, sigma]
p0 = [Rs_est, 1.5e-6, 0.80, Rp_est, 15.0]

bounds = (
    [1.0,   1e-9, 0.3, 10.0, 0.01],   # Lower bounds
    [500.0, 1e-3, 1.0, 1e7,  500000.0]  # Upper bounds
)

popt_E, pct_E, chi2_E, chi2_red_E = fit_model(Z_model_E, p0, bounds, omega, Z_exp)
print("Model E Results:")
columns = ["Rs", "Q1", "n1", "R1", "sigma", "chi2", "chi2_red"]
row_values = [popt_E[0], popt_E[1], popt_E[2], popt_E[3], popt_E[4], chi2_E, chi2_red_E]
row_errors = [pct_E[0],  pct_E[1],  pct_E[2],  pct_E[3],  pct_E[4],  None,    None]
table_E = pd.DataFrame([row_values, row_errors], columns=columns, index=["Values", "Error %"])
display(table_E)

# Calculate smooth fit for plot
omega_smooth_E = np.logspace(np.log10(omega.min()), np.log10(omega.max()), 500)
Z_fit_smooth_E = Z_model_E(omega_smooth_E, *popt_E)
Z_fit_bode_E = Z_model_E(omega, *popt_E)


# plotting
fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
# Plot experimental data vs. Model fit
ax.scatter(Z_exp.real, -Z_exp.imag, color='#1F77B4', edgecolor='k', s=65,
           alpha=0.85, zorder=3, label='Experimental Data (CO$_2$, 800 °C)'
ax.plot(Z_fit_smooth_E.real, -Z_fit_smooth_E.imag, color='darkblue', lw=2.5,
        zorder=4, label='Model E NCLS Fit')

# Label key frequencies
freq_labels = [100.0, 1.0, 0.1, 0.01, 0.001]
for fl in freq_labels:
    idx = np.argmin(np.abs(freqs - fl))
    fx, fy = Z_exp[idx].real, -Z_exp[idx].imag
    lbl = f"{fl:.3f} Hz" if fl < 1.0 else f"{fl:.0f} Hz"
    if fl >= 1000.0:
        lbl = f"{fl/1000:.0f} kHz"

    ax.annotate(lbl, xy=(fx, fy), xytext=(12, -12),
                textcoords='offset points', fontsize=9.5, fontweight='bold',
                color='#333333', arrowprops=dict(arrowstyle='->', color='#555555', lw=0.8))
    ax.scatter(fx, fy, color='gold', edgecolor='black', s=90, zorder=5)
# Styling
ax.set_xlabel(r"$Z'$ — Real Impedance  $[\Omega\cdot\mathrm{cm}^2]$", fontsize=13, fontweight='bold', labelpad=8)
ax.set_ylabel(r"$-Z''$ — Imaginary Impedance  $[\Omega\cdot\mathrm{cm}^2]$", fontsize=13, fontweight='bold', labelpad=8)
ax.set_title("Nyquist Plot & Equivalent Circuit Fit \nFe20Cr25NiNb Alloy in CO$_2$ at 800 °C (48h - at OCP) ",
             fontsize=14, fontweight='bold', pad=15)
ax.set_aspect('equal', adjustable='box')
ax.set_xlim(left=0, right=max(Zr) * 1.28)
ax.set_ylim(bottom=0, top=max(Zi) * 1.28)
ax.grid(True, which='both', linestyle='--', linewidth=0.6, alpha=0.5)
ax.legend(fontsize=11, frameon=True, edgecolor='#cccccc', facecolor='#ffffff', loc='upper left')
plt.tight_layout()

# bode plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 9.5), sharex=True, dpi=300)
mag_exp = np.abs(Z_exp)
phase_exp = data["Phase (deg)"]
mag_fit_E = np.abs(Z_fit_bode_E)
phase_fit_E = np.angle(Z_fit_bode_E, deg=True)

# Subplot 1: Magnitude
ax1.loglog(freqs, mag_exp, 'o', color='#1F77B4', markeredgecolor='k', markersize=7, alpha=0.85, label='Experimental')

ax1.loglog(freqs, mag_fit_E, '-', color='darkblue', lw=2.5, label='Model E Fit')
ax1.set_ylabel(r"$|Z|$ Magnitude  $[\Omega\cdot\mathrm{cm}^2]$", fontsize=12, fontweight='bold')
ax1.set_title("Bode Plot: Magnitude & Phase Spectrum", fontsize=14, fontweight='bold', pad=15)
ax1.grid(True, which="both", linestyle='--', linewidth=0.5, alpha=0.5)
ax1.legend(fontsize=11, loc='lower left')

# Subplot 2: Phase Angle
ax2.semilogx(freqs, phase_exp, 'o', color='#1F77B3', markeredgecolor='k', markersize=7, alpha=0.85, label='Experimental')
ax2.semilogx(freqs, phase_fit_E, '-', color='darkblue', lw=2.5, label='Model E Fit')
ax2.set_xlabel("Frequency  $[Hz]$", fontsize=12, fontweight='bold')
ax2.set_ylabel("Phase Angle  $[deg]$", fontsize=12, fontweight='bold')
ax2.grid(True, which="both", linestyle='--', linewidth=0.5, alpha=0.5)
ax2.legend(fontsize=11, loc='lower right')
plt.tight_layout()




















