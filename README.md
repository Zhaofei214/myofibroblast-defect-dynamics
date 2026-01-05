# Data Archiving

**Author:** Zhaofei Zheng  
**Date:** January 2026

## Introduction
This repository contains MATLAB-based analysis code for paper 
**“Myofibroblasts slow down defect recombination dynamics in mixed cell monolayers.”**

---

## Software Requirements

### MATLAB
The analysis code was tested using **MATLAB R2023b Update 3 (24.2.0.2806996)** on macOS (Sequoia 15.6.1).

### Third-Party Tools and External Software
The analysis relies on the following third-party tools:

- **PIVlab** (MATLAB toolbox)  
  Used for particle image velocimetry analysis of cell velocity fields.  
  https://pivlab.blogspot.com/

- **Cellpose** (external segmentation software)  
  Used for cell segmentation. Cellpose-generated outputs are imported into MATLAB for downstream analysis.  
  https://www.cellpose.org/

- **ImageJ/Fiji** with the **TrackMate** plugin  
  Used for nucleus detection and tracking. TrackMate outputs are exported and analyzed in MATLAB.  
  https://imagej.net/software/fiji/

Third-party software is not included in this repository and must be installed separately by the user.

---

The pipeline can be executed step-by-step or orchestrated using a single entry
point (`run_all.m`).

---

## Data Archiving and Processing Pipeline

### 1. Director Field Calculation and Defect Detection

The main script for director field computation and topological defect detection is:

- `cziprocessing_defectcounting.m`

This script processes raw microscopy data (e.g. `.czi` files), computes local cell
orientations, and identifies topological defects while conserving topological charge.

It calls or depends on the following functions:
- `directorField.m` — computes the director field from cell orientation or texture data.
- `visualizeCellOrientation.m` — visualizes the local cell orientation as a nematic director field.
- `get_all_defects.m` — identifies and locates all topological defects.

**Outputs:**
- `director_data_*.mat`
- `defectData.mat`
- `processed_results.csv`

---

### 2. Defect Type Counting

Defect statistics are analyzed using:

- `countdefecttype.m`

This script counts defects by type or charge (e.g. +1/2 vs −1/2, integer vs half-integer)
as a function of time or experimental condition.

**Inputs:**
- `defectData.mat`

**Outputs:**
- Summary tables and figures saved under `Results/`

---

### 3. Defect Velocity Analysis

Defect motion is quantified using:

- `defectvelocity.m`

This script tracks defect trajectories over time and computes defect velocities
and displacement statistics.

**Inputs:**
- `defectData.mat`
- (Optional) `director_data_*.mat`

**Outputs:**
- Defect velocity statistics and plots under `Results/`

---

### 4. P-Distribution / Angular Statistics

Orientational statistics are computed using:

- `Pdistribution.m`

This script calculates angular or order-parameter distributions (e.g. \(P(\theta)\))
to characterize nematic alignment and orientational order.

**Inputs:**
- `director_data_*.mat`

**Outputs:**
- Distribution plots and summary statistics under `Results/`

---

### 5. Defect Density Decay Analysis

Temporal decay of +1/2 defect density is analyzed using:

- `defect_density_decay.m`

This script reads preprocessed CSV files containing defect counts and time indices,
normalizes defect density by the physical imaging area, and fits the decay using an
exponential model:

\[
\rho(t) = A \exp(-t / B) + C
\]

**Outputs:**
- Publication-quality plots of defect density versus time
- Text files summarizing fitted parameters, confidence intervals, and goodness-of-fit

---

### 6. Color-Coded Theta Maps

Spatial maps of nematic orientation are generated using:

- `color_angle_plot.m`

This script computes the local orientation angle
\(\theta = \mathrm{mod}(\arctan(n_y/n_x), \pi)\)
from the director field and exports color-coded orientation maps with physical units.

**Inputs:**
- `director_data_*.mat`

**Outputs:**
- Theta maps saved under `Results/theta_maps/`

---

### 7. Velocity, Director, and Defect Overlay

Cell velocity fields are combined with director fields and defect locations using:

- `overlay_velocity_director.m`

This script:
- Loads PIV velocity data (`PIVlab.mat`)
- Computes velocity magnitude (µm/h)
- Resizes velocity maps to full-resolution images
- Overlays director fields and defect positions

**Inputs:**
- `PIVlab.mat`
- `director_data_*.mat`
- `defectData.mat`

**Outputs:**
- Velocity magnitude `.mat` files
- High-resolution TIFF overlays in `Results/`

---

### 8. Cell Ratio Verification (Magenta/Green Channels)

Relative cell population ratios are estimated using:

- `cell_ratio_verification.m`

This script thresholds red and green fluorescence channels and computes relative
cell area fractions to verify experimental mixing ratios.

**Inputs:**
- RGB fluorescence images (`.tif` / `.png`)

**Outputs:**
- Cell area ratios printed to the MATLAB command window

---

### 9. Details for using cellPose

Cell segmentation is performed using Cellpose in a Python environment and segmentation outputs are subsequently imported into MATLAB for director and defect analysis. Here is an example of bash command used in this project

conda activate /Users/<your user name>/miniconda3/envs/cellpose
cd ~/Desktop/Cellpose/750_70_03
cellpose --dir . --pretrained_model cyto2 --chan 0 --chan2 0 --save_tif --verbose


### 10. Density, Director, and Defect Overlay

Cell density maps are combined with nematic director fields and topological defect locations using:

- `overlay_density_directors_defects.m`

This script:
- Loads green and red cell density maps (cells/µm²)
- Constructs single-channel (green/red) and combined RGB density backgrounds
- Resamples and overlays nematic director fields
- Overlays ±1/2 topological defect positions
- Exports publication-quality figures (300 dpi)

**Inputs:**
- `segmentation/density_maps_frame*.mat`
- `defectimages/director_data_*.mat`
- `defectData.mat`

**Outputs:**
- `overlay_green_directors_frame*.tif`
- `overlay_red_directors_frame*.tif`
- `overlay_combined_directors_frame*.tif`  
  (saved in `Results/`)

### 11. Velocity Spatial Correlation Analysis

Spatial correlations of cell velocity fields are quantified using:

- `velocity_spatial_correlation_two_conditions.m`

This script computes the isotropic spatial velocity correlation function \(C_v(r)\)
from PIV-derived velocity fields and compares two experimental conditions
(e.g., different cell densities).

**Inputs:**
- `PIVlab.mat` (from two condition folders, e.g. `750_50/` and `500_50/`)
  - `u_filtered`
  - `v_filtered`

**Outputs:**
- Spatial velocity correlation plot (`velocity_correlation.png`) saved under `Results/`
- Printed mean velocity magnitudes for each condition

### 12 Detection on defect and cell dominance

Defects are classified by their local cellular composition (Green vs Pink) and
topological charge (+1/2 or −1/2) using **binarized density maps**.

- `defect_color_dominance.m`

This script:
- Loads defect positions and charges from `defectData.mat`
- Loads per-frame Green/Red density maps and binarizes them (z-score > 0)
- Classifies each defect based on local G/R dominance within a circular neighborhood
- Counts Green(+), Green(−), Pink(+), Pink(−), and Unknown defects per frame
- Produces a clean visualization with four-sided frames and no axis numbers

**Inputs:**
- `defectData.mat`
- `segmentation/density_maps_frame*.mat`
  - `G_cells_per_um2`
  - `R_cells_per_um2`
- `defectimages/director_data_*.mat`

**Outputs:**
- Defect count table (`defect_type_counts.csv`) saved under `Results/`


### 13. Run-All Entry Point

The entire pipeline can be orchestrated using:

- `run_all.m`

This script initializes the repository environment and provides toggle switches
to execute individual analysis steps in sequence.

To run the pipeline:

```matlab
run_all
