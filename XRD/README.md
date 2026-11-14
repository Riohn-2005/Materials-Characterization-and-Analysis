# X-Ray Diffraction (XRD) Pattern Visualization Module

Automated multi-duration X-ray diffraction (XRD) pattern processing and vertical waterfall plotting for tracking crystalline oxide phase evolution during high-temperature isothermal oxidation.

## Project Objective

During high-temperature exposure of austenitic stainless steels in $\text{CO}_2$ environments, surface oxide scales evolve over time from transient iron-rich oxides to protective chromia ($\alpha\text{-Cr}_2\text{O}_3$) and spinel phases ($\text{MnCr}_2\text{O}_4 / \text{FeCr}_2\text{O}_4$). This module automates the alignment and vertical offset stacking of diffractograms across multiple exposure durations ($12\text{ h}$, $48\text{ h}$, and $120\text{ h}$) to visualize relative peak growth, phase emergence, and substrate peak attenuation.

## Technical Equations & Modeling

### 1. Bragg's Law (Crystalline Phase Identification)

Constructive interference from crystal lattice planes of interplanar spacing $d_{hkl}$ occurs at diffraction angles $2\theta$ governed by Bragg's Law:

$$n\lambda = 2 d_{hkl} \sin\theta$$

where:
- $n$ is the diffraction order ($n = 1$)
- $\lambda$ is the incident X-ray wavelength (e.g., $\text{Cu } K\alpha = 1.5406\text{ \AA}$)
- $d_{hkl}$ is the lattice spacing corresponding to Miller indices $(hkl)$
- $\theta$ is the Bragg angle (in degrees)

### 2. Automated Vertical Stacking Offset

To compare multi-duration diffraction patterns without peak overlap while preserving relative intensity scaling, each successive diffractogram $I_k(2\theta)$ is shifted vertically by a constant multiple of the global maximum intensity:

$$I_{\text{shifted}, k}(2\theta) = I_k(2\theta) + k \cdot \Delta I_{\text{offset}}, \qquad \Delta I_{\text{offset}} = 0.5 \times \max_{k, \theta}\left(I_k(2\theta)\right)$$

where $k \in \{0, 1, 2\}$ indexes the oxidation exposure intervals ($12\text{ h}$, $48\text{ h}$, $120\text{ h}$).

## Key Code Implementation

### 1. Data Stacking and Visualizing Code

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load multi-duration diffractogram dataset
df = pd.read_csv('xrd_oxidation_data.csv')

plt.figure(figsize=(10, 8), dpi=150)

# Dynamic vertical stacking offset based on global peak intensity
shift_step = df[['12h', '48h', '120h']].max().max() * 0.5

plt.plot(df['2 Theta_12'], df['12h'], label='12 h', color='black', linewidth=1.2)
plt.plot(df['2 Theta_48'], df['48h'] + shift_step, label='48 h', color='blue', linewidth=1.2)
plt.plot(df['2 Theta_120'], df['120h'] + (2 * shift_step), label='120 h', color='red', linewidth=1.2)

plt.xlabel(r'2$\theta$ (degrees)', fontsize=12)
plt.ylabel('Intensity (a.u.)', fontsize=12)
plt.title(r'XRD Patterns of Austenitic SS (Isothermal Oxidation at 800 °C in CO$_2$)', fontsize=13, pad=12)

plt.grid(True, linestyle='--', alpha=0.5)
plt.minorticks_on()
plt.legend(frameon=True, facecolor='white', framealpha=0.9, loc='upper right')
plt.tight_layout()
plt.show()
```

## Output

<img width="1485" height="1185" alt="image" src="https://github.com/user-attachments/assets/79226a28-5a0a-4443-bbac-17547be0cfc5" />

### 2. Peak Annotation Code

```python



# 1. Load data
df = pd.read_csv('xrd_oxidation_data.csv')

# 2. Set  publication style parameters
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True

fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

# 3. Vertical stacking offset (adjust according to dataset intensity range)
shift_step = df[['12h', '48h', '120h']].max().max() * 0.5

# 4. Plot stacked datasets
ax.plot(df['2 Theta(12h)'], df['12h'], label='12 h', color='black', lw=1.1)
ax.plot(df['2 Theta(48h)'], df['48h'] + shift_step, label='48 h', color='#1f77b4', lw=1.1)
ax.plot(df['2 Theta(120h)'], df['120h'] + (2 * shift_step), label='120 h', color='#d62728', lw=1.1)

# 5. Define Peak Annotations: (2theta_deg, y_peak_coordinate, 'Phase_Label')
# Example placeholder coordinates calibrated to the top shifted curve
annotations = [
    (25.0, 1200, 'Phase A'),
    (35.0, 1450, 'Phase B'),
    (43.5, 1800, 'Substrate'),
    (51.0, 1650, 'Substrate'),
    (65.0, 1300, 'Phase B'),
    (74.5, 1550, 'Substrate'),
]

# 6. Add vertical markers and rotated text labels above peaks
for x, y_pos, label in annotations:
    ax.annotate(
        label,
        xy=(x, y_pos),
        xytext=(x, y_pos + 150),
        ha='center',
        va='bottom',
        fontsize=8,
        rotation=90,
        arrowprops=dict(arrowstyle='->', lw=0.6, color='black')
    )

# 7. Axes formatting, minor ticks, and phase key box
ax.set_xlim(20, 80)
ax.set_xlabel(r'2$\theta$ (degrees)', fontsize=11)
ax.set_ylabel('Intensity (a.u.)', fontsize=11)
ax.set_title('Stacked XRD Patterns Across Isothermal Exposure Durations', fontsize=11, fontweight='bold', pad=12)

ax.minorticks_on()
ax.grid(True, which='major', linestyle='--', alpha=0.4)
ax.legend(frameon=True, facecolor='white', framealpha=0.9, fontsize=9, loc='upper left')

phase_key = 'Substrate: Matrix Peak\nPhase A: Primary Oxide\nPhase B: Secondary Phase'
ax.text(0.98, 0.96, phase_key, transform=ax.transAxes, fontsize=8,
        va='top', ha='right',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))

plt.tight_layout()
plt.show()
```

## Data Availability Note: 
Raw experimental diffraction scans and proprietary alloy phase assignments are excluded from this public repository. To adapt the script to your own dataset, update the annotations list with your indexed $(2\theta, I)$ coordinates and corresponding JCPDS/ICDD phase labels.
