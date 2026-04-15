%% ---- 1. File paths (to change) ----
v1_img_path  = 'YOUR_PATH\V1.nii.gz';
wm_mask_path = 'YOUR_PATH\WM_mask.nii.gz';
out_path     = 'Resultat_theta.nii.gz';

%% ---- 2. Load the images ----
v1_nii = load_nii(v1_img_path);
v1_data = double(v1_nii.img);

wm_nii = load_nii(wm_mask_path);
wm_mask = double(wm_nii.img);

fprintf('  V1      : %s\n', mat2str(size(v1_data)));
fprintf('  WM mask : %s\n', mat2str(size(wm_mask)));

assert(isequal(size(v1_data,1), size(wm_mask,1)) && ...
       isequal(size(v1_data,2), size(wm_mask,2)) && ...
       isequal(size(v1_data,3), size(wm_mask,3)), ...
       'Error with the dimensions');

assert(size(v1_data,4) == 3, ...
    'V1 must have 3 components in the 4th dimension');

%% ---- 3. Compute theta from V1 ----
numerator = sum(v1_data(:,:,:,1:2).^2, 4);   % Vx^2 + Vy^2
denominator = sum(v1_data.^2, 4);            % Vx^2 + Vy^2 + Vz^2

sin_theta_squared = numerator ./ denominator;
theta = asin(sqrt(sin_theta_squared)) * (180 / pi);

%% ---- 4. Apply WM mask ----
wm_binary = double(wm_mask > 0);
theta = theta .* wm_binary;

%% ---- 5. Remove NaN/Inf ----
theta(isnan(theta) | isinf(theta)) = 0;

%% ---- 6. Display stats inside WM only ----
theta_wm = theta(wm_binary > 0);

fprintf('  Theta min  : %.4f deg\n', min(theta_wm));
fprintf('  Theta max  : %.4f deg\n', max(theta_wm));
fprintf('  Theta mean : %.4f deg\n', mean(theta_wm));

%% ---- 7. Save output ----
out_nii = v1_nii;
out_nii.img = single(theta);

fprintf('Saving in: %s\n', out_path);
save_nii(out_nii, out_path);
