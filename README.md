# ANISO-QSM: An Anisotropic Nonlinear Inversion for Susceptibility-separation Optimization algorithm.

## Overview
Inspired by other susceptibility-separation method, the ANISO-QSM algorithm estimates positive and negative susceptibility. 

## Requirements:
- MATLAB
- Toolboxes:
  
	[FANSI toolbox](https://gitlab.com/cmilovic/FANSI-toolbox)

	[Nifti toolbox](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)

	[Image processing toolbox](https://www.mathworks.com/products/image-processing.html)

	[Parallel computing toolbox](https://www.mathworks.com/products/parallel-computing.html)

	[Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html)

## Input:
- Local field map - 3D (Hz)
- GRE magnitude - 4D
- QSM map - 3D (ppm)
- Brain mask - 3D
- R2' map (=R2*-R2) - 3D (s)
- Reference mask - 3D
- Theta map - 3D (degres)

## Preprocess
In order to obtain the input files presented above, several maps must be generated from MRI data. This section describes how each required input is obtained.

### 1. Local Field Map (LocalField.nii)
The local field map is obtained from the phase images through two main preprocessing steps:
- Phase unwrapping (using ROMEO)
- Background field removal (using V-SHARP)
  
These steps can be performed using the SEPIA toolbox (https://sepia-documentation.readthedocs.io/en/latest/getting_started/Installation.html).


### 2. QSM Map (QSM.nii)   
The QSM map is computed from the local field using dipole inversion.
Recommended method: MEDI (Morphology Enabled Dipole Inversion).


### 3. Brain Mask (maskBrain.nii)   
The brain mask is used to restrict computations to brain tissue.
Recommended tool: FSL - use the BET (Brain Extraction Tool) module


### 4. R2′ Map (R2_Prime.nii)  
The R2′ map is defined as: R2′=R2∗−R2
Where :
- R2* = 1/T2* -> obtained from GRE (Gradient Echo) data
- R2 = 1/T2   -> obtained from SE (Spin Echo) data
Both maps are estimated using mono-exponential decay fit 

### 5. Reference Mask (ref_mask.nii)  
The reference mask is used as a reference region for susceptibility normalization.


### 6. Theta Map (theta.nii)
The theta map is computed from DWI (Diffusion Weighted Imaging) data.
A script to compute this map is provided, "tetha_calculation.m".

## How to use:
- Add the required toolboxes to your path
- Add all input data and place them in the same folder and add it to the path
- Update the "baseDir" variable in the script to point to this folder
- Run example_ANISO_QSM.m

## Output:
- X-positive (χ+): positive susceptibility
- X-negative (χ−): negative susceptibility
