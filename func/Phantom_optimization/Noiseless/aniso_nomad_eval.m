function [obj, rmse_p, rmse_m] = aniso_nomad_eval(mu1, mu2, maxOuterIter, saveMinusIter)
% Evaluate ANISO-QSM objective for NOMAD without saving maps.
% Objective = RMSE( ROI-mean X+ ) + RMSE( ROI-mean |X-| )
%
% Uses persistent cache so GT/seg/mask/etc. are loaded once.

persistent cached
if isempty(cached)
    % ------------------- PATHS (YOUR PATHS) -------------------
    cached.fld_path   = '/Users/danielridani/Desktop/Susceptibility_separation/with_FA/Sepia_localfield.nii.gz';
    cached.mag_path   = '/Users/danielridani/Desktop/Susceptibility_separation/with_FA/Magnitude_data.nii.gz';
    cached.msk_path   = '/Users/danielridani/Desktop/Phd_Polytechnique/Project1/osfstorage/Masks/BrainMask.nii.gz';
    cached.R2s_path   = '/Users/danielridani/Desktop/Susceptibility_separation/with_FA/R2_prime.nii.gz';
    cached.csf_path   = '/Users/danielridani/Desktop/Phd_Polytechnique/Project1/osfstorage/Masks/BrainMask.nii.gz';
    cached.theta_path = '/Users/danielridani/Desktop/Phd_Polytechnique/Project1/osfstorage/Noise/theta.nii.gz';

    cached.GTXP_path  = '/Users/danielridani/Desktop/Phd_Polytechnique/Project1/osfstorage/Susceptibility_Separation_Results/Chi_positive.nii.gz';
    cached.GTXN_path  = '/Users/danielridani/Desktop/Susceptibility_separation/with_FA/Chi_negative_with_anisotropy.nii.gz';
    cached.segP_path  = '/Users/danielridani/Desktop/Phd_Polytechnique/Project1/osfstorage/Masks/SegmentedModel.nii.gz';       % labels 1..9
    cached.segN_path  = '/Users/danielridani/Desktop/Phd_Polytechnique/Project1/osfstorage/Masks/white_matter_mask.nii.gz';    % labels 1..10 (must be label map)

    % ------------------- LOAD FIXED DATA ONCE -------------------
    fld_nii   = load_nii(cached.fld_path);
    mag_nii   = load_nii(cached.mag_path);
    msk_nii   = load_nii(cached.msk_path);
    R2s_nii   = load_nii(cached.R2s_path);
    csf_nii   = load_nii(cached.csf_path);
    theta_nii = load_nii(cached.theta_path);

    cached.localField = double(fld_nii.img);
    cached.magn       = double(mag_nii.img);
    cached.mask       = double(msk_nii.img);
    cached.R2s        = double(R2s_nii.img);
    cached.csf        = double(csf_nii.img);
    cached.theta      = double(theta_nii.img);

    GTXP_nii = load_nii(cached.GTXP_path);
    GTXN_nii = load_nii(cached.GTXN_path);
    segP_nii = load_nii(cached.segP_path);
    segN_nii = load_nii(cached.segN_path);

    cached.GTXP = double(GTXP_nii.img);
    cached.GTXN = double(GTXN_nii.img); % negative GT
    cached.segP = double(segP_nii.img);
    cached.segN = double(segN_nii.img);

    % ------------------- CLEAN -------------------
    cached.localField(~isfinite(cached.localField)) = 0;
    cached.magn(~isfinite(cached.magn))             = 0;
    cached.mask(~isfinite(cached.mask))             = 0;
    cached.R2s(~isfinite(cached.R2s))               = 0;
    cached.theta(~isfinite(cached.theta))           = 0;
    cached.theta(isnan(cached.theta))               = 0;

    cached.magn = cached.magn .* cached.mask;

    % convert to Hz if rad/s
    cached.localField = cached.localField / (2*pi);

    % ------------------- Kernel + weights (fixed) -------------------
    cached.vox = [0.64 0.64 0.64];
    cached.kernel = dipole_kernel_fansi(size(cached.mask), cached.vox, 0);

    % magnitude can be 3D or 4D
    if ndims(cached.magn) >= 4
        tmp = sqrt(mean(cached.magn.^2, 4));
    else
        tmp = abs(cached.magn);
    end
    wmap = (tmp ./ max(tmp(:))) .* cached.mask;
    cached.wmap = wmap / max(wmap(:));

    Gm    = gradient_calc(cached.wmap, 1);
    SNR   = 1000;
    noise = 1 / SNR;
    Gm    = max(Gm, noise);
    cached.regweight = min(mean(Gm(:)) ./ Gm, 100);

    % labels
    cached.labelsP = 1:9;
    cached.labelsN = 1:10;

    % GT X- magnitude
    cached.GTXN_pos = -cached.GTXN;

    % fixed scalars (your current ones)
    cached.Drp = 323;
    cached.Drn = 401;

    % fixed hyperparams you want fixed
    cached.alpha1 = 0.01;   % FIXED FOR NOW
    cached.tolUpdate = 0.1; % FIXED
    cached.isGPU = false;
    cached.isPrecond = true;
end

% ------------------- BUILD PARAMS -------------------
params = struct();
params.input         = cached.localField;
params.K             = cached.kernel;
params.weight        = cached.wmap;
params.regweight     = cached.regweight;

params.alpha1        = cached.alpha1;
params.mu1           = mu1;
params.mu2           = mu2;

params.maxOuterIter  = maxOuterIter;
params.saveMinusIter = saveMinusIter;
params.tolUpdate     = cached.tolUpdate;

params.isGPU         = cached.isGPU;
params.isPrecond     = cached.isPrecond;

params.R2s           = cached.R2s;
params.Drp           = cached.Drp;
params.Drn           = cached.Drn;
params.csf           = cached.csf;
params.theta         = cached.theta;

% ------------------- RUN -------------------
out = ANISO_QSM(params);

xp = double(out.x_plus) .* cached.mask;
xm = double(out.x_minus) .* cached.mask; % negative 
xm_pos = abs(xm);

% ------------------- ROI-MEAN RMSE (X+) -------------------
gt_means_p  = [];
rec_means_p = [];
per_roi_rmse_p = [];

for L = cached.labelsP
    idx = (cached.segP == L) & (cached.mask > 0);
    if ~any(idx(:)), continue; end

    gt_mean  = mean(cached.GTXP(idx));
    rec_mean = mean(xp(idx));

    gt_means_p(end+1)  = gt_mean;   
    rec_means_p(end+1) = rec_mean;  

    per_roi_rmse_p(end+1) = sqrt((rec_mean - gt_mean)^2);  % = abs(rec_mean-gt_mean)
end

if isempty(per_roi_rmse_p)
    rmse_p = 1e9;
else
    rmse_p = mean(per_roi_rmse_p);
end


% ------------------- ROI-MEAN "manual RMSE" (X-) in magnitude space -------------------
gt_means_m  = [];
rec_means_m = [];
per_roi_rmse_m = [];

for L = cached.labelsN
    idx = (cached.segN == L) & (cached.mask > 0);
    if ~any(idx(:)), continue; end

    gt_mean  = mean(cached.GTXN_pos(idx));  % GT magnitude
    rec_mean = mean(xm_pos(idx));           % predicted magnitude

    gt_means_m(end+1)  = gt_mean;   
    rec_means_m(end+1) = rec_mean;  

    per_roi_rmse_m(end+1) = sqrt((rec_mean - gt_mean)^2); 
end

if isempty(per_roi_rmse_m)
    rmse_m = 1e9;
else
    rmse_m = mean(per_roi_rmse_m);
end
obj = (rmse_m*0.7+rmse_p*0.3)/2;
end