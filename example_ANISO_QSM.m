%% ========================================================================
%  ANISO-QSM Susceptibility Separation
% ========================================================================

clear
close all
clc

%% ========================================================================
%  1. Load input images
% ========================================================================

fld_nii = load_untouch_nii('path to local field');
mag_nii = load_untouch_nii('path to magnitude data');
msk_nii = load_untouch_nii('path to mask');
R2p_nii = load_untouch_nii('path to R2 prime');
theta_nii = load_untouch_nii('path to theta');
FA_nii    = load_untouch_nii('path to FA');
QSM_nii = load_untouch_nii('path to QSM');


%% ========================================================================
%  2. Convert images to double precision
% ========================================================================

localField = double(fld_nii.img);   % rad/s
magn = double(mag_nii.img);

slope  = double(mag_nii.hdr.dime.scl_slope);
offset = double(mag_nii.hdr.dime.scl_inter);

if slope ~= 0
    magn = magn .* slope + offset;
end
mask       = double(msk_nii.img);
R2p        = double(R2p_nii.img);
theta      = double(theta_nii.img);

FA  = double(FA_nii.img);
QSM = double(QSM_nii.img);


%% ========================================================================
%  3. Remove NaN and Inf values
% ========================================================================

localField(~isfinite(localField)) = 0;
magn(~isfinite(magn))             = 0;
mask(~isfinite(mask))             = 0;
R2p(~isfinite(R2p))               = 0;
theta(~isfinite(theta))           = 0;

FA(~isfinite(FA))   = 0;
QSM(~isfinite(QSM)) = 0;

theta(isnan(theta)) = 0;

magn = magn .* mask;


%% ========================================================================
%  4. Convert local field from rad/s to Hz
% ========================================================================

localField = localField / (2*pi);


%% ========================================================================
%  5. Create dipole kernel and spatial weights
% ========================================================================

vox    = double(fld_nii.hdr.dime.pixdim(2:4)); 
B0_dir = [0 0 1];
B0     = 3;             % Tesla
gyro   = 42.58;         % MHz/T

kernel = dipole_kernel(size(mask), vox);

tmp  = sqrt(mean(magn.^2, 4));
wmap = (tmp ./ max(tmp(:))) .* mask;

Gm    = gradient(wmap);
SNR   = 1000;
noise = 1 / SNR;

Gm = max(Gm, noise);

regweight = min(mean(Gm(:)) ./ Gm, 100);


%% ========================================================================
%  6. Prepare ANISO-QSM parameters
% ========================================================================

params = struct();

% Data terms
params.input  = localField;
params.K      = kernel;
params.weight = wmap / max(wmap(:));

% Regularization and ADMM parameters
params.alpha1        = 0.020801;
params.mu1           = 2080.1;
params.mu2           = 0.5; 
params.maxOuterIter  = 15; 
params.tolUpdate     = 0.1; 
params.regweight     = regweight;

% Solver options
params.isGPU     = false;
params.isPrecond = true;

% Relaxation and susceptibility parameters
params.R2p = R2p;
params.Drp = 140;
params.Drn = 323;

% Masks and orientation maps
params.mask  = mask;
params.theta = theta;

% Initialization maps
params.FA  = FA;
params.QSM = QSM;


%% ========================================================================
%  7. Run ANISO-QSM
% ========================================================================

out = ANISO_QSM(params);


%% ========================================================================
%  8. Prepare output maps
% ========================================================================

chi_ppm_p = single(out.x_plus .* mask);
chi_ppm_m = single(out.x_minus .* mask * -1);
chi_ppm_t = single(chi_ppm_p - chi_ppm_m);


%% ========================================================================
%  9. Save output maps
% ========================================================================

fld_nii.hdr.dime.datatype = 16;
fld_nii.hdr.dime.bitpix   = 32;

fld_nii.img = chi_ppm_p;
save_untouch_nii(fld_nii, 'Chip_ANISO.nii.gz');

fld_nii.img = chi_ppm_m;
save_untouch_nii(fld_nii, 'Chim_ANISO.nii.gz');

fld_nii.img = chi_ppm_t;
save_untouch_nii(fld_nii, 'Chitot_ANISO.nii.gz');

disp('ANISO-QSM finished successfully.')
