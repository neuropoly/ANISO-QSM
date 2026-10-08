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
- [FANSI toolbox](https://gitlab.com/cmilovic/FANSI-toolbox)
- [NIfTI toolbox](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
- [Image Processing Toolbox](https://www.mathworks.com/products/image-processing.html)
- [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html)
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

These steps can be performed using the [SEPIA toolbox](https://sepia-documentation.readthedocs.io/en/latest/getting_started/Installation.html).

### 2. QSM map (`QSM.nii`)
The QSM map is computed from the local field using dipole inversion.

Recommended method:
- **MEDI** (Morphology Enabled Dipole Inversion)

### 3. Brain mask (`maskBrain.nii`)
The brain mask is used to restrict computations to brain tissue.

Recommended tool:
- **FSL BET** (Brain Extraction Tool)

### 4. R2′ map (`R2_Prime.nii`)
The R2′ map is defined as:

`R2′ = R2* - R2`

where:
- `R2* = 1 / T2*`, obtained from GRE data
- `R2 = 1 / T2`, obtained from spin-echo data

Both maps are estimated using a mono-exponential decay fit.

### 5. Reference mask (`ref_mask.nii`)
This mask can be either the whole brain mask or just the CSF mask.

### 6. Theta map (`theta.nii`)
The theta map is computed from diffusion-weighted imaging (DWI) data.

A script is provided:
- `theta_calculation.m`

### 7. Initial Xp and Xn maps (`Xp.nii`, `Xn.nii`)
These maps are obtained from an anisotropy-naive susceptibility-separation algorithm and are used as initialization for ANISO-QSM.

### 8. FA map (`FA.nii`)
The FA map is computed from DWI data and is used as a weighting factor in the anisotropy modeling.

---

## Optimization

To optimize the ANISO-QSM parameters, we recommend using the realistic digital phantom provided in the **Susceptibility-Separation-Phantom** repository together with **NOMAD**, a derivative-free blackbox optimization framework.

### Recommended resources
- [Phantom](https://github.com/neuropoly/Susceptibility-Separation-Phantom)
- [NOMAD](https://www.gerad.ca/en/software/nomad/)

### Installation
Install NOMAD before running the optimization pipeline.  
We recommend using the **Python-based workflow**, since the optimization scripts provided for ANISO-QSM are written in Python and call the MATLAB evaluation functions.

### Optimization modes
Two optimization cases are available:
- **Without noise**
- **With noise**

These allow parameter tuning under ideal simulated conditions or under more realistic noisy conditions.

### Files to modify
Before running the optimization, update the paths in the following files so they match your local installation:

- `aniso_nomad_eval.m`
- `aniso_nomad_save_best.m`
- `nomad_opt_aniso_qsm.py`

### Running the optimization
Once the paths are updated, run:

```bash
python nomad_opt_aniso_qsm.py
```
---
## Usage

1. Add the required toolboxes to the MATLAB path.
2. Place all input files in the same folder.
3. Update the `baseDir` variable in `example_ANISO_QSM.m` so that it points to this folder.
4. Run:

```matlab
example_ANISO_QSM
```
---
## Outputs

ANISO-QSM generates the following output maps:
- χ⁺: positive susceptibility map in ppm
- χ⁻: negative susceptibility map in ppm
