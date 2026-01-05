## bfmatlab (Bio-Formats MATLAB Toolbox)

This directory contains the **Bio-Formats MATLAB toolbox (`bfmatlab`)**, which provides
MATLAB functions for reading and writing a wide range of biological microscopy file
formats (e.g. `.czi`, `.lif`, `.nd2`, `.ome.tiff`).

The toolbox is developed and maintained by the Open Microscopy Environment (OME) and is
used in this repository for loading raw microscopy data prior to downstream processing
(e.g. segmentation, PIV analysis, director field extraction).

### Typical usage in this repository

The toolbox is used to:
- Read raw multi-dimensional microscopy files (Z, T, C)
- Access image metadata (pixel size, time intervals, channels)
- Convert proprietary formats into MATLAB-compatible arrays

Example functions commonly used:
- `bfopen`
- `bfGetReader`
- `bfGetPlane`

### Installation / Setup

This toolbox is **not written by the authors of this repository** and is included for
convenience.

To enable it in MATLAB, add the folder to your path:

```matlab
addpath(genpath('path/to/bfmatlab'));