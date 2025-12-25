# Data Archiving

**Author:** Zhaofei Zheng  
**Date:** November 2025

## Introduction
This repository contains MATLAB-based analysis code for the PNAS manuscript  
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

## Data Archiving and Processing Pipeline

### 1. Director Field Calculation and Defect Detection
The main script for director field computation and defect analysis is:

- `cziprocessing_defectcounting.m`

This script calls the following functions:
- `directorField.m` — computes the director field from cell orientation data.
- `visualizeCellOrientation.m` — visualizes the local cell orientation as a director field.
- `get_all_defects.m` — identifies and locates all topological defects.

---

### 2. Cell Density and Defects
Cell density is further analyzed using outputs from **PIVlab**.  
The file `PIVlab.mat` serves as the input for the MATLAB script:

- `overlayVelocityDensity.m`

This script overlays the cell velocity field with the density map and generates violin plots for comparative analysis.

---

### 3. Cellpose / Segmentation
Cell segmentation is performed using **Cellpose** in a Python environment.  
The following commands are used:

```bash
conda activate /Users/zhaofei/miniconda3/envs/cellpose
cd ~/Desktop/Cellpose/750_70_03
cellpose --dir . --pretrained_model cyto2 --chan 0 --chan2 0 --save_tif --verbose