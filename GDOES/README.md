# GDOES Depth Profile Visualization Module

Glow Discharge Optical Emission Spectroscopy (GDOES) depth profile visualization and analysis for oxide scale thickness and metal/oxide interface characterization. GDOES is a rapid, destructive depth-profiling technique used to map elemental distributions from the surface into the bulk substrate.

## Project Objective

This module provides automated tools for visualizing and analyzing GDOES depth profiles, showing elemental composition as a function of sputtering depth. It supports quantitative estimation of oxide layer thickness, interface width, and elemental depletion/enrichment gradients in high-temperature oxidation and protective coating studies.

## Technical Equations & Modeling

### Depth Calibration

GDOES sputters material layer-by-layer using an argon plasma discharge. Sputtering time is converted to depth via:

$$\text{Depth} = r_{\text{sputter}} \times t$$

where:
- $\text{Depth}$ is the sputtered depth ($\text{nm}$ or $\mu\text{m}$)
- $r_{\text{sputter}}$ is the material-specific calibrated sputter rate ($\text{nm/s}$ or $\mu\text{m/s}$)
- $t$ is the sputtering time ($\text{s}$)

### Concentration Normalization

For quantitative compositional depth profiling, elemental signals are normalized to $100\%$:

$$C_i' = \frac{C_i}{\sum_j C_j} \times 100$$

where $C_i$ is the calibrated concentration of element $i$.

## Visualizer Preview

> **Data Availability Note:** This visualizer was originally developed and validated using high-temperature oxidation datasets of austenitic stainless steel across multiple exposure durations. Because the raw experimental datasets are proprietary research data, they are excluded from this public repository; representative interface screenshots are shown below.

Single Sample depth Profile
<img width="1680" height="1050" alt="Screenshot 2026-11-14 at 1 50 08 PM" src="https://github.com/user-attachments/assets/f6b79cf3-c169-4277-84f4-fb314e723733" />

Time Series Evolution of Cr
<img width="1680" height="1050" alt="Screenshot 2026-11-14 at 1 50 51 PM" src="https://github.com/user-attachments/assets/4a9e3b52-c1a6-48bc-82f4-9a792e9ef959" />

## Expected Input Data Format

The script accepts comma-separated (`.csv`) files containing a depth or time axis alongside elemental concentrations (in $\text{at.}\%$ or $\text{wt.}\%$):

- `Depth (nm)` or `Time (s)`: Sputtered depth or elapsed discharge time
- Element columns: e.g., `Fe`, `Cr`, `Ni`, `O`, `Mn`, `Nb`

```csv
Depth (nm),Fe,Cr,Ni,O
0.0,40.2,30.1,5.1,24.6
20.2,40.4,30.0,5.1,24.5

