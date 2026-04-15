% =========================================================
%  CALCUL DE THETA A PARTIR DE V1
%  Entrées :
%    - V1_nii.gz
%    - SegmentedModel.nii
%  Sortie :
%    - Resultat_theta.nii.gz
% =========================================================


%% ---- 1. File path (To be changed)----
v1_img_path = 'YOUR_PATH\V1.nii.gz';
seg_path    = 'YOUR_PATH\SegmentedModel.nii';
out_path    = 'Resultat_theta.nii.gz';

%% ---- 2. Load the image ----
v1_info  = niftiinfo(v1_img_path);
v1_data  = double(niftiread(v1_info));     % 256 x 320 x 320 x 3

seg_info = niftiinfo(seg_path);
seg_data = double(niftiread(seg_info));    % 256 x 320 x 320

fprintf('  V1  : %s | type : %s\n', mat2str(v1_info.ImageSize),  v1_info.Datatype);
fprintf('  Seg : %s | type : %s\n', mat2str(seg_info.ImageSize), seg_info.Datatype);

assert(isequal(v1_info.ImageSize(1:3), seg_info.ImageSize(1:3)), ...
    'Error with the dimensions');

assert(v1_info.ImageSize(4) == 3, ...
    'V1 must have 3 dimensions');

numerator         = sum(v1_data(:,:,:,1:2) .^ 2, 4);   % Vx² + Vy²
denominator = sum(v1_data .^ 2, 4);          % (Vx²+Vy²+Vz²)²
sin_theta_squared = numerator ./ denominator;

theta = asin(sqrt(sin_theta_squared)) * (180 / pi);

seg_binary = double(seg_data == 8);  
theta = theta .* seg_binary;

theta(isnan(theta) | isinf(theta)) = 0;

fprintf('  Theta min : %.4f deg\n', min(theta(seg_data > 0)));
fprintf('  Theta max : %.4f deg\n', max(theta(seg_data > 0)));
fprintf('  Theta mean : %.4f deg\n', mean(theta(seg_data > 0)));

%% ---- 8. Header ----
out_info                  = v1_info;
out_info.Datatype         = 'single';
out_info.BitsPerPixel     = 32;
out_info.ImageSize        = size(theta);                   
out_info.PixelDimensions  = v1_info.PixelDimensions(1:3);  
%% ---- 9. Save ----
fprintf('Saving in: %s\n', out_path);
niftiwrite(single(theta), out_path, out_info, 'Compressed', true);
