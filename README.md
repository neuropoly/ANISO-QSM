# ANISO-QSM: An Anisotropic Nonlinear Inversion for Susceptibility-separation Optimization algorithm.

## Requirements:
- MATLAB
- [FANSI toolbox](https://gitlab.com/cmilovic/FANSI-toolbox)
- [Nifti toolbox](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image) 

## How to use:
- Add the FANSI toolbox to your path
- Run example_ANISO_QSM.m

## Input:
- Local field map ( after applying phase unwrapping + background field removal) - 3D
- GRE magnitude - 4D
- QSM map - 3D
- Brain mask - 3D
- R2' map (=R2*-R2) - 3D
- CSF mask (as a reference) - 3D
- theta map (from DTI) - 3D

## Output:
- X-positive
- X-negative

