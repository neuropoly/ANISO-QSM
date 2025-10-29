
    %% ========================================================================
    %%  FANSI – Susceptibility Separation (phase-based, 3 T data)
    %% ========================================================================
    addpath(genpath(pwd)); % FANSI toolbox root
    
    % -------- 1. Load your images ------------------------------------------------
    fld_nii = load_nii('/Users/danielridani/Desktop/Phd_Polytechnique/Data_UNF/chi_002/nifti/QSM/Sepia_localfield.nii.gz');
    mag_nii = load_nii('/Users/danielridani/Desktop/Phd_Polytechnique/Data_UNF/chi_002/nifti/Magnitude_data.nii.gz');
    msk_nii = load_nii('/Users/danielridani/Desktop/Phd_Polytechnique/Data_UNF/chi_002/nifti/mask2.nii.gz');
    R2s_nii = load_nii('/Users/danielridani/Desktop/Phd_Polytechnique/Data_UNF/chi_002/nifti/R2_prime_denoised.nii.gz');
    csf_nii = load_nii('/Users/danielridani/Desktop/Phd_Polytechnique/Data_UNF/chi_002/nifti/mask2.nii.gz');
    theta_nii = load_nii('/Users/danielridani/Desktop/Phd_Polytechnique/Data_UNF/chi_002/nifti/registered_theta.nii.gz');                % %% === NEW
    
    localField = double(fld_nii.img); % rad · s⁻¹ (ω)
    magn = double(mag_nii.img); % echo-1 magnitude
    mask = double(msk_nii.img);
    R2s = double(R2s_nii.img);
    csf = double(csf_nii.img);
    theta = double(theta_nii.img);                                        

    localField(~isfinite(localField)) = 0;
    magn(~isfinite(magn)) = 0;
    mask(~isfinite(mask)) = 0;
    R2s(~isfinite(R2s)) = 0;
    theta(~isfinite(theta)) = 0; 
    theta(~isnan(theta)) = 0; 

    magn=magn.*mask;

    % -------- 3. Convert unit -------------------------
    localField = localField / (2*pi); 

    % -------- 4. Dipole kernel ---------------------------------------------------
    vox = [0.7 0.7 0.7]; % mm
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
    params.K = kernel;
    params.weight = wmap / max(wmap(:));
    params.alpha1 = 0.1;
    params.mu1 = 10;
    params.mu2 = 1;
    params.isGPU = false;
    params.isPrecond = true;
    params.maxOuterIter = 5; 
    params.saveMinusIter = 5;
    params.tolUpdate = 0.1;
    params.regweight = regweight;
    params.R2s = R2s;
    params.Drp = 120;
    params.Drn = 120;
    params.csf = csf;
    params.theta = theta;    
    % -------- 6. Run susceptibility separation ------------------------------------
    out = ANISO_QSM(params);

    % ----------------  SAVE  --------------------
    chi_ppm_p = single(out.x_plus  .* mask);
    chi_ppm_m = single(out.x_minus .* mask*-1);

    fld_nii.hdr.dime.datatype = 16; fld_nii.hdr.dime.bitpix = 32;
    fld_nii.img = chi_ppm_p; save_nii(fld_nii, 'Chip_ANISO_chi_002.nii.gz' );
    fld_nii.img = chi_ppm_m; save_nii(fld_nii,'Chim_ANISO_chi_002.nii.gz' );
      