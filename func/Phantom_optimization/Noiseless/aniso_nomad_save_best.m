function aniso_nomad_save_best(mu1, mu2, maxOuterIter, saveMinusIter, out_dir)
% Run ANISO-QSM once with best params and save x_plus + x_minus.

% Paths (same as eval)
fld_path   = '/Users/danielridani/Desktop/Susceptibility_separation/with_FA/Sepia_localfield.nii.gz';
mag_path   = '/Users/danielridani/Desktop/Susceptibility_separation/with_FA/Magnitude_data.nii.gz';
msk_path   = '/Users/danielridani/Desktop/Phd_Polytechnique/Project1/osfstorage/Masks/BrainMask.nii.gz';
R2s_path   = '/Users/danielridani/Desktop/Susceptibility_separation/with_FA/R2_prime.nii.gz';
csf_path   = '/Users/danielridani/Desktop/Phd_Polytechnique/Project1/osfstorage/Masks/BrainMask.nii.gz';
theta_path = '/Users/danielridani/Desktop/Phd_Polytechnique/Project1/osfstorage/Noise/theta.nii.gz';

fld_nii   = load_nii(fld_path);
mag_nii   = load_nii(mag_path);
msk_nii   = load_nii(msk_path);
R2s_nii   = load_nii(R2s_path);
csf_nii   = load_nii(csf_path);
theta_nii = load_nii(theta_path);

localField = double(fld_nii.img);
magn       = double(mag_nii.img);
mask       = double(msk_nii.img);
R2s        = double(R2s_nii.img);
csf        = double(csf_nii.img);
theta      = double(theta_nii.img);

localField(~isfinite(localField)) = 0;
magn(~isfinite(magn))             = 0;
mask(~isfinite(mask))             = 0;
R2s(~isfinite(R2s))               = 0;
theta(~isfinite(theta))           = 0;
theta(isnan(theta))               = 0;

magn = magn .* mask;
localField = localField / (2*pi);

vox    = [0.64 0.64 0.64];
kernel = dipole_kernel_fansi(size(mask), vox, 0);

if ndims(magn) >= 4
    tmp = sqrt(mean(magn.^2, 4));
else
    tmp = abs(magn);
end
wmap = (tmp ./ max(tmp(:))) .* mask;

Gm    = gradient_calc(wmap, 1);
SNR   = 1000;
noise = 1 / SNR;
Gm    = max(Gm, noise);
regweight = min(mean(Gm(:)) ./ Gm, 100);

params = struct();
params.input         = localField;
params.K             = kernel;
params.weight        = (wmap / max(wmap(:)));
params.regweight     = regweight;

params.alpha1        = 0.01;   % FIXED FOR NOW
params.mu1           = mu1;
params.mu2           = mu2;

params.maxOuterIter  = maxOuterIter;
params.saveMinusIter = saveMinusIter;
params.tolUpdate     = 0.1;

params.isGPU         = false;
params.isPrecond     = true;

params.R2s           = R2s;
params.Drp           = 323;
params.Drn           = 401;
params.csf           = csf;
params.theta         = theta;

out = ANISO_QSM(params);

xp = double(out.x_plus) .* mask;
xm = double(out.x_minus) .* mask;

if ~exist(out_dir, 'dir'), mkdir(out_dir); end

ref_nii = msk_nii;

ref_nii.img = xp;
save_nii(ref_nii, fullfile(out_dir, 'x_plus_best.nii.gz'));

ref_nii.img = xm;
save_nii(ref_nii, fullfile(out_dir, 'x_minus_best.nii.gz'));
end