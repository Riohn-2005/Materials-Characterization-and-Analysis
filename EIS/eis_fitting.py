#!/usr/bin/env python3
"""
EIS Equivalent Circuit Fitting Module

This module provides tools for fitting electrochemical impedance spectroscopy (EIS) data
to various equivalent circuit models using Constant Phase Elements (CPE) and Warburg elements.

Physical Background:
- Impedance Z = Z' - jZ'' where Z' is real part, Z'' is imaginary part
- CPE Impedance: Z_CPE = 1 / [Q(jω)^n], where Q is CPE constant, n is dispersion coefficient (0 ≤ n ≤ 1)
- Warburg Impedance (semi-infinite diffusion): Z_W = σ / √(jω)
- Finite Warburg Impedance: Z_W = R_W * tanh[(jωτ_W)^n_W] / (jωτ_W)^n_W

Author: Materials Characterization Toolkit
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from typing import Tuple, Dict, Callable, Optional


# ============================================================================
# EQUIVALENT CIRCUIT MODELS
# ============================================================================

def z_model_cpe_parallel(omega: np.ndarray, Rs: float, Q1: float, n1: float, R1: float) -> complex:
    """
    Simplest model: Rs - (CPE1 // R1)
    
    Circuit: Solution resistance Rs in series with parallel CPE1 and R1
    
    Z_total = Rs + Z_parallel
    Z_parallel = 1 / (1/Z_CPE + 1/R1) = 1 / [Q(jω)^n + 1/R1]
    
    Parameters:
        omega: Angular frequency (rad/s)
        Rs: Solution resistance (Ohm·cm²)
        Q1: CPE constant (s^n / Ohm·cm²)
        n1: CPE dispersion coefficient (0 ≤ n1 ≤ 1)
        R1: Polarization resistance (Ohm·cm²)
    
    Returns:
        Complex impedance (Ohm·cm²)
    """
    j = 1j
    Z_parallel = 1.0 / (Q1 * (j * omega)**n1 + 1.0 / R1)
    return Rs + Z_parallel


def z_model_warburg_parallel(omega: np.ndarray, Rs: float, Q1: float, n1: float, R1: float, sigma: float) -> complex:
    """
    Model E: Rs - CPE1 // (R1 - W)
    
    Circuit: Rs in series with parallel CPE1 and branch containing R1 and Warburg element in series
    
    Z_W = σ / √(jω)  (semi-infinite Warburg diffusion)
    
    Parameters:
        omega: Angular frequency (rad/s)
        Rs: Solution resistance (Ohm·cm²)
        Q1: CPE constant (s^n / Ohm·cm²)
        n1: CPE dispersion coefficient
        R1: Charge transfer resistance (Ohm·cm²)
        sigma: Warburg coefficient (Ohm·s^-0.5·cm²)
    
    Returns:
        Complex impedance (Ohm·cm²)
    """
    j = 1j
    Z_W = sigma / np.sqrt(j * omega)
    Z_branch = R1 + Z_W
    Z_parallel = 1.0 / (Q1 * (j * omega)**n1 + 1.0 / Z_branch)
    return Rs + Z_parallel


def z_model_series_parallel(omega: np.ndarray, Rs: float, Q1: float, n1: float, R1: float, 
                            Q2: float, n2: float, R2: float) -> complex:
    """
    Model F: Rs - (CPE1 // R1) - (CPE2 // R2)
    
    Circuit: Two parallel R-CPE branches in series with Rs
    
    Parameters:
        omega: Angular frequency (rad/s)
        Rs: Solution resistance (Ohm·cm²)
        Q1, Q2: CPE constants (s^n / Ohm·cm²)
        n1, n2: CPE dispersion coefficients
        R1, R2: Polarization resistances (Ohm·cm²)
    
    Returns:
        Complex impedance (Ohm·cm²)
    """
    j = 1j
    Z_p1 = 1.0 / (Q1 * (j * omega)**n1 + 1.0 / R1)
    Z_p2 = 1.0 / (Q2 * (j * omega)**n2 + 1.0 / R2)
    return Rs + Z_p1 + Z_p2


def z_model_nested_parallel(omega: np.ndarray, Rs: float, Q1: float, n1: float, R1: float,
                            Q2: float, n2: float, R2: float) -> complex:
    """
    Model J: Rs - (CPE1 // (R1 + (CPE2 // R2)))
    
    Circuit: Nested parallel configuration
    
    Parameters:
        omega: Angular frequency (rad/s)
        Rs: Solution resistance (Ohm·cm²)
        Q1, Q2: CPE constants
        n1, n2: CPE dispersion coefficients
        R1, R2: Resistances (Ohm·cm²)
    
    Returns:
        Complex impedance (Ohm·cm²)
    """
    j = 1j
    Z_inner = 1.0 / (Q2 * (j * omega)**n2 + 1.0 / R2)
    Z_branch = R1 + Z_inner
    Z_outer = 1.0 / (Q1 * (j * omega)**n1 + 1.0 / Z_branch)
    return Rs + Z_outer


def z_model_nested_diffusion(omega: np.ndarray, Rs: float, Q1: float, n1: float, R1: float,
                             Q2: float, n2: float, R2: float, sigma: float) -> complex:
    """
    Model B: Nested Inner Boundary Diffusion (Rs - CPE1 // (R1 - CPE2 // (R2 - W)))
    
    Circuit: Nested structure with Warburg diffusion in innermost branch
    
    Parameters:
        omega: Angular frequency (rad/s)
        Rs: Solution resistance (Ohm·cm²)
        Q1, Q2: CPE constants
        n1, n2: CPE dispersion coefficients
        R1, R2: Resistances (Ohm·cm²)
        sigma: Warburg coefficient
    
    Returns:
        Complex impedance (Ohm·cm²)
    """
    j = 1j
    Z_W = sigma / np.sqrt(j * omega)
    Z_inner = 1.0 / (Q2 * (j * omega)**n2 + 1.0 / (R2 + Z_W))
    Z_outer = 1.0 / (Q1 * (j * omega)**n1 + 1.0 / (R1 + Z_inner))
    return Rs + Z_outer


def z_model_finite_warburg(omega: np.ndarray, Rs: float, Q1: float, n1: float, R1: float,
                           Q2: float, n2: float, R2: float, Rw: float, tauW: float, nw: float = 0.5) -> complex:
    """
    Model D: Nested Inner Bounded Finite Diffusion
    
    Uses finite-length (short-circuit) Warburg impedance:
    Z_W = R_W * tanh[(jωτ_W)^n_W] / (jωτ_W)^n_W
    
    Parameters:
        omega: Angular frequency (rad/s)
        Rs: Solution resistance (Ohm·cm²)
        Q1, Q2: CPE constants
        n1, n2: CPE dispersion coefficients
        R1, R2: Resistances (Ohm·cm²)
        Rw: Warburg resistance (Ohm·cm²)
        tauW: Warburg time constant (s)
        nw: Warburg exponent (typically 0.5)
    
    Returns:
        Complex impedance (Ohm·cm²)
    """
    j = 1j
    s_w = (j * omega * tauW) ** nw
    Z_W = Rw * np.tanh(s_w) / s_w
    Z_inner = 1.0 / (Q2 * (j * omega)**n2 + 1.0 / (R2 + Z_W))
    Z_outer = 1.0 / (Q1 * (j * omega)**n1 + 1.0 / (R1 + Z_inner))
    return Rs + Z_outer


def z_model_triple_parallel(omega: np.ndarray, Rs: float, Q1: float, n1: float, R1: float,
                            Q2: float, n2: float, R2: float, Q3: float, n3: float, R3: float) -> complex:
    """
    Model C: Rs - CPE1 // (R1 - (R2 // CPE2) - (R3 // CPE3))
    
    Circuit: Outer CPE1 parallel to series R1 with two parallel R-CPE blocks
    
    Parameters:
        omega: Angular frequency (rad/s)
        Rs: Solution resistance (Ohm·cm²)
        Q1, Q2, Q3: CPE constants
        n1, n2, n3: CPE dispersion coefficients
        R1, R2, R3: Resistances (Ohm·cm²)
    
    Returns:
        Complex impedance (Ohm·cm²)
    """
    j = 1j
    Z_block3 = 1.0 / (Q3 * (j * omega)**n3 + 1.0 / R3)
    Z_block2 = 1.0 / (Q2 * (j * omega)**n2 + 1.0 / R2)
    Z_lower_branch = R1 + Z_block2 + Z_block3
    Z_outer_parallel = 1.0 / (Q1 * (j * omega)**n1 + 1.0 / Z_lower_branch)
    return Rs + Z_outer_parallel


def z_model_triple_series(omega: np.ndarray, Rs: float, Q1: float, n1: float, R1: float,
                          Q2: float, n2: float, R2: float, Q3: float, n3: float, R3: float) -> complex:
    """
    Model I: Series Voigt Circuit (Rs + (CPE1 // R1) + (CPE2 // R2) + (CPE3 // R3))
    
    Circuit: Three parallel R-CPE blocks in series with Rs
    
    Parameters:
        omega: Angular frequency (rad/s)
        Rs: Solution resistance (Ohm·cm²)
        Q1, Q2, Q3: CPE constants
        n1, n2, n3: CPE dispersion coefficients
        R1, R2, R3: Resistances (Ohm·cm²)
    
    Returns:
        Complex impedance (Ohm·cm²)
    """
    j = 1j
    Z_p1 = 1.0 / (Q1 * (j * omega)**n1 + 1.0 / R1)
    Z_p2 = 1.0 / (Q2 * (j * omega)**n2 + 1.0 / R2)
    Z_p3 = 1.0 / (Q3 * (j * omega)**n3 + 1.0 / R3)
    return Rs + Z_p1 + Z_p2 + Z_p3


# ============================================================================
# FITTING FUNCTIONS
# ============================================================================

def fit_eis_model(model_func: Callable, p0: list, bounds: Tuple[list, list],
                 omega: np.ndarray, Z_exp: np.ndarray,
                 lin_indices: Optional[set] = None) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """
    Fit an EIS model to experimental impedance data using weighted least squares.
    
    Physical weighting: residuals are weighted by |Z| to properly account for
    magnitude variations across frequency spectrum.
    
    Parameters:
        model_func: Impedance model function
        p0: Initial parameter guesses
        bounds: Tuple of (lower_bounds, upper_bounds)
        omega: Angular frequency array (rad/s)
        Z_exp: Experimental complex impedance array (Ohm·cm²)
        lin_indices: Set of parameter indices to keep in linear space (not log-transformed)
                     Typically CPE exponents (n values)
    
    Returns:
        popt: Optimized parameters
        pct_err: Percentage errors for each parameter
        chi2: Chi-squared goodness of fit
        chi2_red: Reduced chi-squared
    """
    if lin_indices is None:
        lin_indices = set()
    
    n_params = len(p0)
    
    def residuals(opt_params):
        """Weighted residuals function for optimization"""
        params = [np.exp(val) if i not in lin_indices else val for i, val in enumerate(opt_params)]
        Zm = model_func(omega, *params)
        w = np.abs(Z_exp)  # Weight by impedance magnitude
        return np.concatenate([(Z_exp.real - Zm.real) / w, (Z_exp.imag - Zm.imag) / w])
    
    # Convert to log-space for positive parameters (except exponents)
    opt_p0 = [val if i in lin_indices else np.log(val) for i, val in enumerate(p0)]
    opt_bounds = (
        [b if i in lin_indices else np.log(b) for i, b in enumerate(bounds[0])],
        [b if i in lin_indices else np.log(b) for i, b in enumerate(bounds[1])]
    )
    
    # Multi-start optimization to escape local minima
    best_res = None
    best_chi2 = np.inf
    
    starts = [opt_p0]
    # Vary CPE exponents in start points
    for n_val in [0.75, 0.85]:
        s = list(opt_p0)
        for idx in lin_indices:
            s[idx] = n_val
        starts.append(s)
    
    for x0 in starts:
        x0 = np.clip(x0, opt_bounds[0], opt_bounds[1])
        try:
            res = least_squares(residuals, x0, bounds=opt_bounds, method='trf',
                              ftol=1e-15, xtol=1e-15, gtol=1e-15)
            chi2_val = np.sum(res.fun**2)
            if chi2_val < best_chi2:
                best_chi2 = chi2_val
                best_res = res
        except Exception:
            continue
    
    # Convert parameters back to linear space
    opt_x = best_res.x
    popt = np.array([opt_x[i] if i in lin_indices else np.exp(opt_x[i]) for i in range(n_params)])
    
    # Calculate statistics and error propagation
    chi2 = np.sum(best_res.fun**2)
    n_points = len(omega)
    dof = 2 * n_points - n_params
    chi2_red = chi2 / dof
    
    # Covariance matrix and standard errors
    J = best_res.jac
    s2 = chi2 / max(dof, 1)
    cov_opt = np.linalg.pinv(J.T @ J) * s2
    std_opt = np.sqrt(np.maximum(0, np.diag(cov_opt)))
    
    # Calculate percentage errors
    perr = np.zeros(n_params)
    for i in range(n_params):
        perr[i] = std_opt[i] if i in lin_indices else popt[i] * std_opt[i]
    
    pct_err = (perr / popt) * 100.0
    
    return popt, pct_err, chi2, chi2_red


# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

def plot_nyquist(Zr: np.ndarray, Zi: np.ndarray, Z_fit: Optional[np.ndarray] = None,
                 freqs: Optional[np.ndarray] = None, title: str = "Nyquist Plot",
                 save_path: Optional[str] = None) -> plt.Figure:
    """
    Plot Nyquist diagram (Z' vs -Z'')
    
    Parameters:
        Zr: Real impedance (Ohm·cm²)
        Zi: Negative imaginary impedance (Ohm·cm²)
        Z_fit: Fitted complex impedance (optional)
        freqs: Frequency array for annotation (optional)
        title: Plot title
        save_path: Path to save figure (optional)
    
    Returns:
        matplotlib Figure object
    """
    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
    
    # Plot experimental data
    ax.plot(Zr, Zi, 'o', label='Experimental Data', markersize=4, alpha=0.7)
    
    # Plot fit if provided
    if Z_fit is not None:
        ax.plot(Z_fit.real, -Z_fit.imag, '-', label='Fit', linewidth=2, color='red')
    
    # Annotate frequency points if provided
    if freqs is not None:
        for i, f in enumerate(freqs[::len(freqs)//10]):  # Annotate ~10 points
            ax.annotate(f'{f:.0e}', xy=(Zr[i*10], Zi[i*10]), fontsize=8, alpha=0.5)
    
    ax.set_xlabel("Z' (Ohm·cm²)", fontsize=12, fontweight='bold')
    ax.set_ylabel("-Z'' (Ohm·cm²)", fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend()
    ax.axis('equal')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


def plot_bode(freqs: np.ndarray, Z_exp: np.ndarray, Z_fit: Optional[np.ndarray] = None,
              title: str = "Bode Plot", save_path: Optional[str] = None) -> plt.Figure:
    """
    Plot Bode diagram (|Z| and phase vs frequency)
    
    Parameters:
        freqs: Frequency array (Hz)
        Z_exp: Experimental complex impedance
        Z_fit: Fitted complex impedance (optional)
        title: Plot title
        save_path: Path to save figure (optional)
    
    Returns:
        matplotlib Figure object
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), dpi=150, sharex=True)
    
    # Magnitude plot
    mag_exp = np.abs(Z_exp)
    ax1.semilogx(freqs, mag_exp, 'o', label='Experimental', markersize=4, alpha=0.7)
    
    if Z_fit is not None:
        mag_fit = np.abs(Z_fit)
        ax1.semilogx(freqs, mag_fit, '-', label='Fit', linewidth=2, color='red')
    
    ax1.set_ylabel('|Z| (Ohm·cm²)', fontsize=12, fontweight='bold')
    ax1.set_title(title, fontsize=13, fontweight='bold')
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.legend()
    
    # Phase plot
    phase_exp = -np.angle(Z_exp) * 180 / np.pi
    ax2.semilogx(freqs, phase_exp, 'o', label='Experimental', markersize=4, alpha=0.7)
    
    if Z_fit is not None:
        phase_fit = -np.angle(Z_fit) * 180 / np.pi
        ax2.semilogx(freqs, phase_fit, '-', label='Fit', linewidth=2, color='red')
    
    ax2.set_xlabel('Frequency (Hz)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Phase (°)', fontsize=12, fontweight='bold')
    ax2.grid(True, linestyle='--', alpha=0.5)
    ax2.legend()
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


# ============================================================================
# DATA LOADING
# ============================================================================

def load_eis_data(filepath: str, cycle: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load EIS data from CSV file.
    
    Expected columns: 'Frequency (Hz)', 'Re(Z) normalized (Ohm.cm²)', '-Im(Z) normalized (Ohm.cm²)'
    Optional column: 'Cycle' for multi-cycle experiments
    
    Parameters:
        filepath: Path to CSV file
        cycle: Cycle number to extract (if multi-cycle data)
    
    Returns:
        freqs: Frequency array (Hz)
        Zr: Real impedance (Ohm·cm²)
        Zi: Negative imaginary impedance (Ohm·cm²)
    """
    df = pd.read_csv(filepath)
    
    if cycle is not None and 'Cycle' in df.columns:
        df = df[df['Cycle'] == cycle]
    
    freqs = df["Frequency (Hz)"].values
    Zr = df["Re(Z) normalized (Ohm.cm²)"].values
    Zi = df["-Im(Z) normalized (Ohm.cm²)"].values
    
    return freqs, Zr, Zi


# ============================================================================
# MAIN EXAMPLE USAGE
# ============================================================================

def main():
    """Example usage of EIS fitting module"""
    
    # Load data (replace with your file path)
    # freqs, Zr, Zi = load_eis_data('example_eis.csv')
    
    # For demonstration, generate synthetic data
    print("Generating synthetic EIS data for demonstration...")
    freqs = np.logspace(0, 6, 50)  # 1 Hz to 1 MHz
    omega = 2 * np.pi * freqs
    
    # Generate synthetic impedance using Model F
    Rs_true, Q1_true, n1_true, R1_true = 10.0, 1e-5, 0.9, 1000.0
    Q2_true, n2_true, R2_true = 1e-4, 0.8, 5000.0
    
    Z_synthetic = z_model_series_parallel(omega, Rs_true, Q1_true, n1_true, R1_true,
                                          Q2_true, n2_true, R2_true)
    
    # Add noise
    np.random.seed(42)
    noise_level = 0.02 * np.abs(Z_synthetic)
    Z_noisy = Z_synthetic + noise_level * (np.random.randn(len(Z_synthetic)) + 
                                          1j * np.random.randn(len(Z_synthetic)))
    
    Zr = Z_noisy.real
    Zi = -Z_noisy.imag
    Z_exp = Zr - 1j * Zi
    
    # Initial parameter estimates
    Rs_est = Zr[0]
    Rp_total = max(Zr[-1] - Rs_est, 10.0)
    
    p0 = [Rs_est, 1e-6, 0.85, Rp_total * 0.3, 1e-5, 0.75, Rp_total * 0.7]
    bounds = (
        [1.0, 1e-9, 0.3, 1.0, 1e-9, 0.3, 1.0],
        [100.0, 1e-2, 1.0, 1e6, 1e-2, 1.0, 1e7]
    )
    
    # Fit model
    print("Fitting Model F (Series-Parallel R-CPE)...")
    popt, pct_err, chi2, chi2_red = fit_eis_model(
        z_model_series_parallel, p0, bounds, omega, Z_exp, lin_indices={2, 5}
    )
    
    # Print results
    print("\n--- Fitting Results ---")
    print(f"Chi-squared: {chi2:.4e}")
    print(f"Reduced Chi-squared: {chi2_red:.4e}")
    print("\nParameters:")
    names = ['Rs', 'Q1', 'n1', 'R1', 'Q2', 'n2', 'R2']
    for name, val, err in zip(names, popt, pct_err):
        print(f"{name:6s}: {val:.4e} ± {err:.2f}%")
    
    # Calculate fitted impedance
    Z_fit = z_model_series_parallel(omega, *popt)
    
    # Plot results
    print("\nGenerating plots...")
    plot_nyquist(Zr, Zi, Z_fit, freqs, title="Nyquist Plot - Model F Fit")
    plot_bode(freqs, Z_exp, Z_fit, title="Bode Plot - Model F Fit")
    
    plt.show()
    
    return popt, chi2_red


if __name__ == "__main__":
    main()
