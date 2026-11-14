# Raman Spectroscopy Peak Deconvolution Module

Automated baseline subtraction, Savitzky-Golay smoothing, and constrained multi-peak Pseudo-Voigt deconvolution for phase identification and semi-quantitative oxide scale analysis in Raman spectroscopy.

## Project Objective

High-temperature oxide scales formed on Cr-Ni austenitic alloys exhibit overlapping vibrational bands from corundum-type oxides ($\alpha\text{-Cr}_2\text{O}_3$, $\alpha\text{-Fe}_2\text{O}_3$) and spinel phases ($\text{Fe}_3\text{O}_4$, $\text{MnCr}_2\text{O}_4$, $\text{FeCr}_2\text{O}_4$). This module isolates fluorescence background signals, deconvolutes overlapping Raman modes within the $150\text{--}1000\text{ cm}^{-1}$ spectral window, and computes integrated peak area fractions for each oxide phase.

## Technical Equations & Modeling

### 1. Asymmetric Least Squares (ALS) Baseline Correction

Fluorescence backgrounds are estimated iteratively by minimizing a penalized least-squares objective function:

$$S = \sum_{i} w_i (y_i - z_i)^2 + \lambda \sum_{i} (\Delta^2 z_i)^2$$

where:
- $y_i$ is the raw measured intensity and $z_i$ is the fitted smooth baseline
- $\lambda$ controls baseline smoothness ($\lambda = 10^6$)
- $w_i$ are asymmetric weights updated iteratively using $w_i = p$ if $y_i > z_i$ and $w_i = 1 - p$ if $y_i \le z_i$ ($p = 0.001$)

### 2. Pseudo-Voigt Line-Shape Profile

Each Raman vibrational mode is modeled as a linear combination of Gaussian ($G$) and Lorentzian ($L$) profiles sharing a common center $\nu_0$ and half-width at half-maximum $\gamma$:

$$I(\nu) = A \left[ \eta \, L(\nu) + (1 - \eta) \, G(\nu) \right]$$

where:

$$G(\nu) = \exp\left(-\ln(2)\left(\frac{\nu - \nu_0}{\gamma}\right)^2\right), \qquad L(\nu) = \frac{1}{1 + \left(\frac{\nu - \nu_0}{\gamma}\right)^2}$$

- $A$ is peak amplitude, $\nu_0$ is peak position ($\text{cm}^{-1}$), $\eta \in [0, 1]$ is the Lorentzian mixing fraction, and $\text{FWHM} = 2\gamma$.


## Key Code Implementation

### Baseline & Pseudo-Voigt Formulation

```python
def baseline_als(y, lam=1e6, p=0.001, niter=15):
    L = len(y)
    D = eye(L, format='csc')
    D = D[1:] - D[:-1]
    D = D[1:] - D[:-1]
    w = np.ones(L)
    for _ in range(niter):
        W = csc_matrix((w, (np.arange(L), np.arange(L))), shape=(L, L))
        Z = W + lam * D.T * D
        z = spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y <= z)
    return z

def pseudo_voigt(x, amp, cen, eta, gamma):
    g = np.exp(-np.log(2) * ((x - cen) / gamma)**2)
    l = 1 / (1 + ((x - cen) / gamma)**2)
    return amp * (eta * l + (1 - eta) * g)
```
```python
initial_guesses = [
    # Low Wavenumber Region: [Amplitude, Center (cm-1), Eta, Gamma]
    15,  307, 0.3, 5,    # P1: α-Fe2O3
    8,   351, 0.3, 5,    # P2: α-Cr2O3 (Weak)
    5,   398, 0.3, 6,    # P3: α-Fe2O3 (Weak)

    # Main Chromia Region
    25,  532, 0.4, 8,    # P4: α-Cr2O3 Shoulder (Eg)
    140, 562, 0.4, 4.5,  # P5: α-Cr2O3 Principal Peak (A1g)

    # Spinel Region
    18,  620, 0.4, 8,    # P6: α-Cr2O3 Eg mode
    30,  656, 0.5, 14,   # P7: Fe3O4 / Iron Spinel
    45,  686, 0.5, 18    # P8: MnCr2O4 / FeCr2O4 Spinel
]
```

## Output

<img width="1089" height="989" alt="image" src="https://github.com/user-attachments/assets/134efb03-07d3-43ba-8497-bdf690935203" />
