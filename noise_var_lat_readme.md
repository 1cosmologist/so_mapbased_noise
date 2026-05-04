# SO LAT Polarization Noise Variance Maps

This document describes the methodology for creating polarized noise variance maps for the Simons Observatory (SO) Large Aperture Telescope (LAT) for both the **wide** and **delensing** scan strategies.

## Overview

The noise variance maps are generated from telescope observation schedules, instrument NET (Noise Equivalent Temperature), and relative hits maps from simulated scan strategies. These maps provide spatially-varying noise estimates for polarization (Q/U) measurements.

Two scan strategies are supported:
- **wide** (`lat_wide_1_el`): Main wide-area survey covering ~40% of sky.
- **delens_wide** (`lat_delens_wide`): Delensing fields overlapping with the wide survey, covering a smaller, deeper patch.

The total telescope-year budget is the same for both strategies (see schedule below). The two surveys can be combined with inverse-variance weighting for a given split fraction using `code/wide-delens_mix.py`.

## Input Data

### Relative Hits Maps

**Wide survey** (`lat_wide_1_el`):
- LF band: `resources/so_lat_lf_relhits_2048.fits`
- MF band: `resources/so_lat_mf_relhits_2048.fits`
- UHF band: `resources/so_lat_uhf_relhits_2048.fits`

**Delensing survey** (`lat_delens_wide`):
- LF band: `resources/so_lat_delens_wide_lf_relhits_2048.fits`
- MF band: `resources/so_lat_delens_wide_mf_relhits_2048.fits`
- UHF band: `resources/so_lat_delens_wide_uhf_relhits_2048.fits`

All hits maps are at nside=2048.

### Instrument Parameters

Parameters are taken from `resources/instr_params_lat-wide_channels.yaml` (wide) and `resources/instr_params_lat-delens_channels.yaml` (delensing).

**Note**: The knee frequency and alpha parameters differ for Intensity (I) and Polarization (P).

| Channel | Frequency (GHz) | Beam FWHM (arcmin) | ℓ_knee (I) | α (I) | ℓ_knee (P) | α (P) |
|---------|-----------------|--------------------|------------|-------|------------|-------|
| LF027   | 27              | 7.4                | 500        | -3.5  | 700        | -1.4  |
| LF039   | 39              | 5.1                | 500        | -3.5  | 700        | -1.4  |
| MF093   | 93              | 2.2                | 2100       | -3.5  | 700        | -1.4  |
| MF145   | 145             | 1.4                | 3000       | -3.5  | 700        | -1.4  |
| HF225   | 225             | 1.0                | 3800       | -3.5  | 700        | -1.4  |
| HF280   | 280             | 0.9                | 3800       | -3.5  | 700        | -1.4  |

### Detector Sensitivity

NET values used in the variance map generation:

| Frequency (GHz) | NET_telescope (μK√s) |
|-----------------|----------------------|
| 27              | 35                   |
| 39              | 18                   |
| 93              | 7.8                  |
| 145             | 8.4                  |
| 225             | 14.1                 |
| 280             | 35.4                 |

These values are assumed to have the detector yield factor folded in.

### Efficiency Factors
- **Observation efficiency**: `obs_eff = 0.30` (Standard LAT efficiency factor, applied identically to both wide and delensing surveys)

---

## Telescope Observation Schedule

### Baseline Schedule (nominal SO 2025 + ASO from mid-2026, through 2035)

The same total tube-year budget is used for both scan strategies. The fraction of time allocated to each strategy per frequency band is controlled externally (see `code/wide-delens_mix.py`).

Number of tubes deployed per year (2025–2035):

| Band  | '25  | '26  | '27  | '28  | '29  | '30  | '31  | '32  | '33  | '34  | '35  | **Total (tube-yrs)** |
|-------|------|------|------|------|------|------|------|------|------|------|------|----------------------|
| f027  | 0    | 1    | 1    | 1    | 1    | 1    | 1    | 1    | 1    | 1    | 1    | **10.0**             |
| f039  | 0    | 1    | 1    | 1    | 1    | 1    | 1    | 1    | 1    | 1    | 1    | **10.0**             |
| f093  | 4    | 8    | 8    | 8    | 8    | 8    | 8    | 8    | 8    | 8    | 8    | **84.0**             |
| f145  | 4    | 8    | 8    | 8    | 8    | 8    | 8    | 8    | 8    | 8    | 8    | **84.0**             |
| f225  | 2    | 4    | 4    | 4    | 4    | 4    | 4    | 4    | 4    | 4    | 4    | **42.0**             |
| f280  | 2    | 4    | 4    | 4    | 4    | 4    | 4    | 4    | 4    | 4    | 4    | **42.0**             |

Each variance map is computed assuming **all** telescope-years are assigned to that survey. When splitting the budget between wide and delensing, variance maps are rescaled by `1/fraction` at run time.

---

## Noise Variance Calculation

### Formula

The polarization noise variance per pixel is calculated as:

$$\sigma^2_{\text{pix}} = \frac{(\sqrt{2} \cdot \text{NET})^2}{t_{\text{obs,pix}} \cdot \eta_{\text{obs}}}$$

Where:
- $\sqrt{2}$ factor converts from NET (temperature) to NEQ/U (polarization)
- $\text{NET}$ is the detector array noise equivalent temperature
- $t_{\text{obs,pix}}$ is the observation time per pixel in seconds
- $\eta_{\text{obs}}$ is the observation efficiency

### Observation Time per Pixel

$$t_{\text{obs,pix}} = \frac{t_{\text{total}} \cdot h_{\text{rel,pix}}}{\sum h_{\text{rel}}}$$

Where:
- $t_{\text{total}}$ is the total observation time (years $\times$ sec/year)
- $h_{\text{rel,pix}}$ is the relative hits at each pixel

---

## Output Variance Maps

### File Naming Convention

**Wide survey:**
```
so_lat_wide_1_el_f{freq}_noise_var_{tel_yrs:.2f}tel-yrs.fits
```

**Delensing survey:**
```
so_lat_delens_wide_f{freq}_noise_var_{tel_yrs:.2f}tel-yrs.fits
```

### Surveys
- **wide** (`lat_wide_1_el`): Main wide-area survey.
- **delens_wide** (`lat_delens_wide`): Delensing fields overlapping with the wide survey.

### Generated Files

#### Wide Survey

| File | Tel-yrs | Description |
|------|---------|-------------|
| `so_lat_wide_1_el_f027_noise_var_10.00tel-yrs.fits` | 10.00 | 27 GHz variance map |
| `so_lat_wide_1_el_f039_noise_var_10.00tel-yrs.fits` | 10.00 | 39 GHz variance map |
| `so_lat_wide_1_el_f093_noise_var_84.00tel-yrs.fits` | 84.00 | 93 GHz variance map |
| `so_lat_wide_1_el_f145_noise_var_84.00tel-yrs.fits` | 84.00 | 145 GHz variance map |
| `so_lat_wide_1_el_f225_noise_var_42.00tel-yrs.fits` | 42.00 | 225 GHz variance map |
| `so_lat_wide_1_el_f280_noise_var_42.00tel-yrs.fits` | 42.00 | 280 GHz variance map |

#### Delensing Survey

| File | Tel-yrs | Description |
|------|---------|-------------|
| `so_lat_delens_wide_f027_noise_var_10.00tel-yrs.fits` | 10.00 | 27 GHz variance map |
| `so_lat_delens_wide_f039_noise_var_10.00tel-yrs.fits` | 10.00 | 39 GHz variance map |
| `so_lat_delens_wide_f093_noise_var_84.00tel-yrs.fits` | 84.00 | 93 GHz variance map |
| `so_lat_delens_wide_f145_noise_var_84.00tel-yrs.fits` | 84.00 | 145 GHz variance map |
| `so_lat_delens_wide_f225_noise_var_42.00tel-yrs.fits` | 42.00 | 225 GHz variance map |
| `so_lat_delens_wide_f280_noise_var_42.00tel-yrs.fits` | 42.00 | 280 GHz variance map |

Each map encodes the white-noise variance per pixel (uK²) assuming the full telescope-year budget is assigned to that survey.

## Polarization Map Depth Summary

Estimated polarization map depths (mean over observed sky, full tel-yr budget):

### Wide Survey

| Channel | Map-depth (μK-arcmin) |
|---------|----------------------|
| f027    | 45                   |
| f039    | 23                   |
| f093    | 3.4                  |
| f145    | 3.7                  |
| f225    | 8.8                  |
| f280    | 22                   |

### Delensing Survey

The delensing field covers a smaller sky area than the wide survey; depths are lower (better) in the overlapping region due to the concentrated scan strategy.

| Channel | Map-depth (μK-arcmin, deepest pixel) |
|---------|--------------------------------------|
| f027    | 17                                  |
| f039    | 8.9                                  |
| f093    | 1.3                                 |
| f145    | 1.4                                 |
| f225    | 3.3                                 |
| f280    | 8.3                                  |

> **Note**: When the budget is split between wide and delensing with fraction `f`, the effective depth at any pixel degrades as `depth / sqrt(fraction)` for single-survey pixels, or improves via ILC combination for pixels covered by both footprints. See `code/wide-delens_mix.py` for the combination logic.

## Map Depth Visualizations

### 27 GHz
![f027 map depth](figures/f27_lat_mapdepth.png)

### 93 GHz
![f093 map depth](figures/f93_lat_mapdepth.png)

### 145 GHz
![f145 map depth](figures/f145_lat_mapdepth.png)

