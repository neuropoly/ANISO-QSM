    %% ========================================================================
    %%  FANSI – Susceptibility Separation
    %% ========================================================================
    clear; clc;
    % -------- 1. Change for the input data file direction ----------------------------------------
    baseDir = 'YOUR_PATH';
    % ---------------------------------------------------------------------------------------------
    
    fld_nii   = load_nii(fullfile(baseDir, 'LocalField.nii'));
    mag_nii   = load_nii(fullfile(baseDir,  'Magnitude_data.nii'));
    msk_nii   = load_nii(fullfile(baseDir, 'maskBrain.nii'));
    R2p_nii   = load_nii(fullfile(baseDir,  'R2_Prime.nii')); 
    ref_nii   = load_nii(fullfile(baseDir,  'ref_mask.nii')); % can be the same as msk_nii
    theta_nii = load_nii(fullfile(baseDir,  'theta.nii'));
    FA_nii    = load_nii(fullfile(baseDir, 'FA.nii'));
    QSM_nii   = load_nii(fullfile(baseDir, 'QSM.nii'));
    Xp_nii    = load_nii(fullfile(baseDir, 'Xp.nii')); % QSM can also be used here
    Xn_nii    = load_nii(fullfile(baseDir, 'Xn.nii')); % QSM can also be used here
   

    QSM = double(QSM_nii.img);
    localField = double(fld_nii.img); % rad · s⁻¹ (ω)
    magn = double(mag_nii.img); 
    mask = double(msk_nii.img);
    R2p = double(R2p_nii.img);
    ref = double(ref_nii.img);
    theta = double(theta_nii.img);  
    Xp = double(Xp_nii.img);
    Xn = double(Xn_nii.img);
    FA = double(FA_nii.img);

    localField(~isfinite(localField)) = 0;
    magn(~isfinite(magn)) = 0;
    mask(~isfinite(mask)) = 0;
    R2p(~isfinite(R2p)) = 0;
    theta(~isfinite(theta)) = 0; 
    theta(isnan(theta)) = 0; 
    Xp(~isfinite(Xp)) = 0;
    Xn(~isfinite(Xn)) = 0;
    FA(~isfinite(FA)) = 0;

    magn=magn.*mask;
    FA = (FA - min(FA(:))) ./ (max(FA(:)) - min(FA(:)));
    FA = FA .* mask;
    % -------- 3. Convert unit -------------------------
    localField = localField / (2*pi); 

    % -------- 4. Dipole kernel ---------------------------------------------------
    vox = [0.64 0.64 0.64]; 
    B0_dir = [0, 0, 1];
    B0 = 3; % Tesla
    gyro = 42.58; % MHz/T
    kernel = dipole_kernel_fansi(size(mask), vox, 0);

    tmp = sqrt(mean(magn.^2, 4));
    wmap = (tmp ./ max(tmp(:))) .* mask;

    Gm = gradient_calc(wmap, 1);
    SNR = 1000;
    noise = 1 / SNR;
    Gm = max(Gm, noise);
    regweight = min(mean(Gm(:)) ./ Gm,100); 
    
    % -------- 5. Prepare the parameter struct ------------------------------
   
    params = struct();
    params.input = localField;
    params.ref = ref;
    params.theta = theta;    
    params.qsm = QSM;
    params.Xp = Xp;
    params.Xn = Xn;
    params.FA = FA;
    params.K = kernel;
    params.weight = wmap / max(wmap(:));
    params.alpha1 = 0.01;
    params.mu1 = 0.1;
    params.mu2 = 1;
    params.isGPU = false;
    params.isPrecond = true;
    params.maxOuterIter = 100;
    params.tolUpdate = 0.1;
    params.regweight = regweight;
    params.R2p = R2p;
    params.Drp = 323;
    params.Drn = 401;
  
 
    % -------- 6. Run susceptibility separation ------------------------------------
    out = ANISO_QSM(params);

   % ----------------  SAVE  --------------------
    chi_ppm_p = single(out.x_plus  .* mask);
    chi_ppm_m = single(out.x_minus .* mask*-1);

    fld_nii.hdr.dime.datatype = 16; fld_nii.hdr.dime.bitpix = 32;

    % -------- 7. Change the image name accordingly ------------------------------------
    fld_nii.img = chi_ppm_p; save_nii(fld_nii, 'Xp.nii.gz' );
    fld_nii.img = chi_ppm_m; save_nii(fld_nii,'Xm.nii.gz' ); 

      
