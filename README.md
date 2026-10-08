# ANISO-QSM

**ANISO-QSM** is an anisotropic nonlinear inversion method for susceptibility separation that estimates **positive susceptibility** (χ⁺) and **negative susceptibility** (χ⁻). **ANISO-QSM** takes into account **susceptibility anisotropy** using information derived from **diffusion-weighted imaging (DWI)** data.

---

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Input Data](#input-data)
- [Preprocessing](#preprocessing)
- [Optimization](#optimization)
- [Usage](#usage)
- [Outputs](#outputs)

---

## Overview

ANISO-QSM is inspired by existing susceptibility-separation methods and extends them by incorporating anisotropy information into the inversion framework.

The algorithm uses the following MRI-derived maps:
- Local field map
- GRE magnitude
- QSM map
- Brain mask
- R2′ map
- Theta map
- FA map

The **Xp** and **Xn** maps are susceptibility-separation results obtained from anisotropy-naive algorithms and are used as an initial starting point for ANISO-QSM. The **FA** map is derived from DWI data and is used as a weighting factor in the anisotropy modeling.

---

## Requirements

### Software
- MATLAB

### MATLAB Toolboxes
- [NIfTI toolbox](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
- [Image Processing Toolbox](https://www.mathworks.com/products/image-processing.html)
- [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html) (optional)
- [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html)

---

## Input Data

The algorithm requires the following input files:

| Input | Description | Dimension | Unit |
|---|---|---:|---:|
| `LocalField.nii` | Local field map | 3D | Hz |
| `Magnitude_data.nii` | GRE magnitude | 4D | - |
| `QSM.nii` | Quantitative susceptibility map | 3D | ppm |
| `maskBrain.nii` | Brain mask | 3D | - |
| `R2_Prime.nii` | R2′ map (`R2* - R2`) | 3D | s<sup>-1</sup> |
| `theta.nii` | Theta map | 3D | degrees |
| `FA.nii` | Fractional anisotropy map derived from DWI | 3D | - |

---

## Preprocessing

Several maps must be generated from MRI data before running ANISO-QSM.

### 1. Local field map (`LocalField.nii`)
The local field map is obtained from phase images using:
- **Phase unwrapping** with **ROMEO**
- **Background field removal** with **V-SHARP**

Recommended tool: [SEPIA toolbox](https://sepia-documentation.readthedocs.io/en/latest/getting_started/Installation.html).

### 2. QSM map (`QSM.nii`)
The QSM map is computed from the local field using dipole inversion.

Recommended method:
- **MEDI** (Morphology Enabled Dipole Inversion)
Recommended tool: [SEPIA toolbox](https://sepia-documentation.readthedocs.io/en/latest/getting_started/Installation.html).

### 3. Brain mask (`maskBrain.nii`)
The brain mask is used to restrict computations to brain tissue.

Recommended tool: **FSL BET** (Brain Extraction Tool)

### 4. R2′ map (`R2_Prime.nii`)
The R2′ map is defined as:

`R2′ = R2* - R2`

where:
- `R2* = 1 / T2*`, obtained from GRE data
- `R2 = 1 / T2`, obtained from spin-echo data

Both maps are estimated using a mono-exponential decay fit.
Recommended tool: [SEPIA toolbox](https://sepia-documentation.readthedocs.io/en/latest/getting_started/Installation.html).

### 5. FA map (`FA.nii`)
The FA map is computed from DWI data and is used as a weighting factor in the anisotropy modeling.

### 6. Theta map (`theta.nii`)
The theta map is computed from diffusion-weighted imaging (DWI) data.


## Usage

1. Add the required toolboxes to the MATLAB path.
2. Update the paths in `example_ANISO_QSM.m` and the desired parameters (the current parameter values are those identified as optimal).
3. Run:

```matlab
example_ANISO_QSM
```
---
## Outputs

ANISO-QSM generates the following output maps:
- χ⁺: positive susceptibility map in ppm
- χ⁻: negative susceptibility map in ppm
- χtotal: total susceptibility map in ppm
