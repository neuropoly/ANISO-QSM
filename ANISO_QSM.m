function out = ANISO_QSM(params, varargin)
totalTime = tic;

%% ========================================================================
%  1. Parse inputs
% ========================================================================

[phase, alpha, N, Kernel, mu, mu2, maxOuterIter, ...
 tolUpdate, regweight, W, ~, isPrecond, R2p, Drp, Drn, ...
 mask, theta, FA, QSM] = parse_inputs_aniso(params, varargin{:});


%% ========================================================================
%  2. Run isotropic optimization
% ========================================================================

fprintf('\n[ANISO-QSM] Computing isotropic susceptibility maps for initialization...\n');

out_iso = run_aniso_solver( ...
    phase, alpha, N, Kernel, mu, mu2, maxOuterIter, ...
    tolUpdate, regweight, W, isPrecond, R2p, Drp, Drn, ...
    mask, theta, FA, QSM, [], [], ...
    'zero', 'isotropic', true);

chi_plus_iso  = out_iso.x_plus;
chi_minus_iso = out_iso.x_minus;

fprintf('[ANISO-QSM] Isotropic initialization completed.\n');


%% ========================================================================
%  3. Run anisotropic optimization
% ========================================================================

fprintf('\n[ANISO-QSM] Running anisotropic susceptibility separation...\n');

out = run_anisotropic_solver( ...
    phase, alpha, N, Kernel, mu, mu2, maxOuterIter, ...
    tolUpdate, regweight, W, isPrecond, R2p, Drp, Drn, ...
    mask, theta, FA, QSM, chi_plus_iso, chi_minus_iso);

fprintf('[ANISO-QSM] Susceptibility separation completed successfully.\n');


%% ========================================================================
%  4. Finalize output
% ========================================================================

out.totalTime = toc(totalTime);

end


%% =========================================================================
%  Anisotropic solver
% =========================================================================

function out = run_anisotropic_solver(phase, alpha, N, Kernel, mu, mu2, ...
                                      maxOuterIter, tolUpdate, ...
                                      regweight, W, isPrecond, R2p, Drp, Drn, ...
                                      mask, theta, FA, QSM, chi_plus_iso, chi_minus_iso)

solverSettings = { ...
    maxOuterIter,  Drn, chi_plus_iso, chi_minus_iso, 'zero',  'anisotropic', true; ...
    5, Drp, [],           [],            'theta', 'isotropic',       false};

solverOutput = cell(size(solverSettings, 1), 1);

for settingIndex = 1:size(solverSettings, 1)

    Dr = solverSettings{settingIndex, 2};

    solverOutput{settingIndex} = run_aniso_solver( ...
        phase, alpha, N, Kernel, mu, mu2, solverSettings{settingIndex, 1}, ...
        tolUpdate, regweight, W, isPrecond, R2p, Dr, Dr, ...
        mask, theta, FA, QSM, solverSettings{settingIndex, 3}, solverSettings{settingIndex, 4}, ...
        solverSettings{settingIndex, 5}, solverSettings{settingIndex, 6}, solverSettings{settingIndex, 7});

end

out = solverOutput{end};

out.x_plus = solverOutput{1}.x_plus;

out.deltaX_xp = solverOutput{1}.deltaX;
out.X0_xp     = solverOutput{1}.X0;

out.deltaX_xn = solverOutput{end}.deltaX;
out.X0_xn     = solverOutput{end}.X0;

out.iter_xp = solverOutput{1}.iter;
out.iter_xn = solverOutput{end}.iter;


end


%% =========================================================================
%  Internal solver
% =========================================================================

function out = run_aniso_solver(phase, alpha, N, Kernel, mu, mu2, maxOuterIter, ...
                                tolUpdate, regweight, W, isPrecond, R2p, Drp, Drn, ...
                                mask, theta, FA, QSM, Xp, Xn, ...
                                c2Mode, initMode, doPrint)

solverTime = tic;


%% ========================================================================
%  1. Prepare anisotropy and scaling terms
% ========================================================================

FA = (FA - min(FA(:))) ./ (max(FA(:)) - min(FA(:)));

switch lower(c2Mode)

    case 'zero'
        c2 = 0;

    case 'theta'
        c2 = cosd(theta).^2;
    otherwise
        error('c2Mode must be either ''zero'' or ''theta''.')

end



Wy = ((W .* phase ./ (W + mu2)));

B0   = 3;
gyro = 42.58;
W2   = B0 * gyro;
mu3 =0.01;

%% ========================================================================
%  2. Initialize ADMM variables
% ========================================================================

z1_dx = zeros(N, 'single');
z1_dy = zeros(N, 'single');
z1_dz = zeros(N, 'single');

s1_dx = zeros(N, 'single');
s1_dy = zeros(N, 'single');
s1_dz = zeros(N, 'single');

z2_dx = zeros(N, 'single');
z2_dy = zeros(N, 'single');
z2_dz = zeros(N, 'single');

s2_dx = zeros(N, 'single');
s2_dy = zeros(N, 'single');
s2_dz = zeros(N, 'single');

z3_dx = zeros(N, 'single');
z3_dy = zeros(N, 'single');
z3_dz = zeros(N, 'single');

s3_dx = zeros(N, 'single');
s3_dy = zeros(N, 'single');
s3_dz = zeros(N, 'single');


%% ========================================================================
%  3. Initialize susceptibility maps
% ========================================================================

switch lower(initMode)

    case 'anisotropic'
        xp = Xp;
        dX = Xn * -1;
        X0 = QSM;

    case 'isotropic'
        xp = QSM;
        dX = QSM;
        X0 = QSM;

    otherwise
        error('initMode must be either ''anisotropic'' or ''isotropic''.')

end


%% ========================================================================
%  4. Initialize data split variable
% ========================================================================

z2 = isPrecond * Wy + (~isPrecond) * zeros(N, 'single');
s2 = zeros(N, 'single');

alpha_over_mu = alpha / mu;


%% ========================================================================
%  5. Fourier gradient operators
% ========================================================================

[k1, k2, k3] = ndgrid(0:N(1)-1, 0:N(2)-1, 0:N(3)-1);

E1 = 1 - exp(2i .* pi .* k1 / N(1));
E2 = 1 - exp(2i .* pi .* k2 / N(2));
E3 = 1 - exp(2i .* pi .* k3 / N(3));

E1t = conj(E1);
E2t = conj(E2);
E3t = conj(E3);

EE2 = E1t .* E1 + E2t .* E2 + E3t .* E3;


%% ========================================================================
%  6. Main ADMM loop
% ========================================================================

if doPrint
    fprintf('%3s\t%10s\n', 'Iter', 'Update');
end

for t = 1:maxOuterIter

     %% --------------------------------------------------------------------
    % 6.1 Update chi-plus
    % ---------------------------------------------------------------------

    tx = E1t .* fftn(z1_dx - s1_dx);
    ty = E2t .* fftn(z1_dy - s1_dy);
    tz = E3t .* fftn(z1_dz - s1_dz);

    xp_prev = xp;

    xp = real(ifftn(( ...
          mu  * (tx + ty + tz) ...
        + mu2 * W2 * conj(Kernel) .* fftn(z2 - s2) ...
        - mu3*mu2 * (W2^2) * abs(Kernel).^2 .* fftn(c2 .* dX + X0) ...
        + fftn(Drp .* (Drn .* (c2 .* dX + X0) + R2p)) ...
        ) ./ ( ...
          eps ...
        + mu  * EE2 ...
        + mu2 * (W2^2) * abs(Kernel).^2 ...
        + Drp .* Drp ...
        )));

    clear tx ty tz


    %% --------------------------------------------------------------------
    % 6.2 Update deltaX
    % ---------------------------------------------------------------------

    tx = E1t .* fftn(z2_dx - s2_dx);
    ty = E2t .* fftn(z2_dy - s2_dy);
    tz = E3t .* fftn(z2_dz - s2_dz);

    dX_prev = dX;

    dX = real(ifftn(( ...
          mu  * (tx + ty + tz) ...
        + mu2 * W2 * conj(Kernel) .* fftn(c2 .* (z2 - s2)) ...
        - mu3*mu2 * (W2^2) * conj(Kernel) .* fftn(c2 .* real(ifftn(Kernel .* fftn(xp + X0)))) ...
        + fftn((Drn .* c2) .* (Drp .* xp - R2p - Drn .* X0)) ...
        ) ./ ( ...
          eps ...
        + mu  * EE2 ...
        + mu2 * (W2^2) * abs(Kernel).^2 ...
        + Drn .* Drn ...
        ))) .* FA;

    clear tx ty tz


    %% --------------------------------------------------------------------
    % 6.3 Update X0
    % ---------------------------------------------------------------------

    tx = E1t .* fftn(z3_dx - s3_dx);
    ty = E2t .* fftn(z3_dy - s3_dy);
    tz = E3t .* fftn(z3_dz - s3_dz);

    X0_prev = X0;

    X0 = real(ifftn(( ...
          mu  * (tx + ty + tz) ...
        + mu2 * W2 * conj(Kernel) .* fftn(z2 - s2) ...
        - mu3*mu2 * (W2^2) * abs(Kernel).^2 .* fftn(xp + c2 .* dX) ...
        + fftn(Drn .* (Drp .* xp - R2p - Drn .* (c2 .* dX))) ...
        ) ./ ( ...
          eps ...
        + mu  * EE2 ...
        + mu2 * (W2^2) * abs(Kernel).^2 ...
        + Drn .* Drn ...
        )));

    clear tx ty tz


    %% --------------------------------------------------------------------
    % 6.4 Check update
    % ---------------------------------------------------------------------

    x_update = 100*norm([xp(:) - xp_prev(:); ...
                           dX(:) - dX_prev(:); ...
                           X0(:) - X0_prev(:)]) ...
                   / (norm([xp(:); dX(:); X0(:)]) + eps);

    if doPrint
        fprintf('%3d\t%10.4f\n', t, x_update);
    end

    if x_update < tolUpdate
        break
    end


    %% --------------------------------------------------------------------
    % 6.5 TV shrinkage for chi-plus
    % ---------------------------------------------------------------------

    Fx = fftn(xp);

    x_dx = real(ifftn(E1 .* Fx));
    x_dy = real(ifftn(E2 .* Fx));
    x_dz = real(ifftn(E3 .* Fx));

    z1_dx = max(abs(x_dx + s1_dx) - regweight(:,:,:,1) * alpha_over_mu, 0) .* sign(x_dx + s1_dx);
    z1_dy = max(abs(x_dy + s1_dy) - regweight(:,:,:,2) * alpha_over_mu, 0) .* sign(x_dy + s1_dy);
    z1_dz = max(abs(x_dz + s1_dz) - regweight(:,:,:,3) * alpha_over_mu, 0) .* sign(x_dz + s1_dz);

    s1_dx = s1_dx + x_dx - z1_dx;
    s1_dy = s1_dy + x_dy - z1_dy;
    s1_dz = s1_dz + x_dz - z1_dz;

    clear Fx x_dx x_dy x_dz


    %% --------------------------------------------------------------------
    % 6.6 TV shrinkage for deltaX
    % ---------------------------------------------------------------------

    Fx = fftn(dX);

    x_dx = real(ifftn(E1 .* Fx));
    x_dy = real(ifftn(E2 .* Fx));
    x_dz = real(ifftn(E3 .* Fx));

    z2_dx = max(abs(x_dx + s2_dx) - regweight(:,:,:,1) * alpha_over_mu, 0) .* sign(x_dx + s2_dx);
    z2_dy = max(abs(x_dy + s2_dy) - regweight(:,:,:,2) * alpha_over_mu, 0) .* sign(x_dy + s2_dy);
    z2_dz = max(abs(x_dz + s2_dz) - regweight(:,:,:,3) * alpha_over_mu, 0) .* sign(x_dz + s2_dz);

    s2_dx = s2_dx + x_dx - z2_dx;
    s2_dy = s2_dy + x_dy - z2_dy;
    s2_dz = s2_dz + x_dz - z2_dz;

    clear Fx x_dx x_dy x_dz


    %% --------------------------------------------------------------------
    % 6.7 TV shrinkage for X0
    % ---------------------------------------------------------------------

    Fx = fftn(X0);

    x_dx = real(ifftn(E1 .* Fx));
    x_dy = real(ifftn(E2 .* Fx));
    x_dz = real(ifftn(E3 .* Fx));

    z3_dx = max(abs(x_dx + s3_dx) - regweight(:,:,:,1) * alpha_over_mu, 0) .* sign(x_dx + s3_dx);
    z3_dy = max(abs(x_dy + s3_dy) - regweight(:,:,:,2) * alpha_over_mu, 0) .* sign(x_dy + s3_dy);
    z3_dz = max(abs(x_dz + s3_dz) - regweight(:,:,:,3) * alpha_over_mu, 0) .* sign(x_dz + s3_dz);

    s3_dx = s3_dx + x_dx - z3_dx;
    s3_dy = s3_dy + x_dy - z3_dy;
    s3_dz = s3_dz + x_dz - z3_dz;

    clear Fx x_dx x_dy x_dz


    %% --------------------------------------------------------------------
    % 6.8 Nonlinear data consistency update
    % ---------------------------------------------------------------------

    chi_minus = c2 .* dX + X0;

    pred_total = real(ifftn(Kernel .* fftn(xp + chi_minus)));

    numerator = W .* sin(z2 - phase) ...
              + mu2 .* z2 ...
              - mu2 .* (pred_total + s2);

    denominator = W .* cos(z2 - phase) ...
                + mu2 ...
                + eps;

    z2 = z2 - numerator ./ denominator;

    s2 = s2 + pred_total - z2;


    %% --------------------------------------------------------------------
    % 6.9 Update Drp and Drn
    % ---------------------------------------------------------------------

    A = [xp(:), -chi_minus(:)];

    coef = (A' * A + 1e-8 * eye(2)) \ (A' * R2p(:));

    Drp = max(coef(1), 0);
    Drn = max(coef(2), 0);

end


%% ========================================================================
%  7. Final post-processing
% ========================================================================

chi_minus = c2 .* dX + X0;

mask    = mask > 0;
sum_map = xp + chi_minus;
beta = mean(sum_map(mask), 'all');
xp        = xp - beta / 2;
chi_minus = chi_minus - beta / 2;

xp(xp < 0) = 0;
chi_minus(chi_minus > 0) = 0;



%% ========================================================================
%  8. Output structure
% ========================================================================

out.time = toc(solverTime);

out.x_plus  = xp;
out.deltaX  = dX;
out.X0      = X0;
out.x_minus = chi_minus;

out.iter       = t;
out.Drp_final  = Drp;
out.Drn_final  = Drn;

end


%% =========================================================================
%  Input parser
% =========================================================================

function [phase, alpha, N, Kernel, mu, mu2, maxOuterIter, tolUpdate, ...
          regweight, W, isGPU, isPrecond, R2p, Drp, Drn, mask, theta, ...
          FA, QSM] = parse_inputs_aniso(params, varargin)

%% -------------------------------------------------------------------------
% Required fields
% --------------------------------------------------------------------------

if ~isstruct(params)
    error('Provide a struct params with required fields.')
end

if isfield(params, 'input')
    phase = params.input;
else
    error('Provide params.input.')
end

if isfield(params, 'alpha1')
    alpha = params.alpha1;
else
    error('Provide params.alpha1.')
end

if isfield(params, 'R2p')
    R2p = params.R2p;
else
    error('Provide params.R2p.')
end

if isfield(params, 'mask')
    mask = params.mask;
else
    error('Provide params.mask.')
end

if isfield(params, 'Drp')
    Drp = params.Drp;
else
    error('Provide params.Drp.')
end

if isfield(params, 'Drn')
    Drn = params.Drn;
else
    error('Provide params.Drn.')
end

if isfield(params, 'theta')
    theta = params.theta;
else
    error('Provide params.theta.')
end

if isfield(params, 'FA')
    FA = params.FA;
else
    error('Provide params.FA.')
end

if isfield(params, 'QSM')
    QSM = params.QSM;
else
    error('Provide params.QSM.')
end

N = size(phase);


%% -------------------------------------------------------------------------
% Default values
% --------------------------------------------------------------------------

defaultMu        = 100 * alpha;
defaultMu2       = 1.0;
defaultNoOuter   = 150;
defaultTol       = 0.1;
defaultRegweight = ones([N 3]);
defaultW         = ones(N);


%% -------------------------------------------------------------------------
% Input parser
% --------------------------------------------------------------------------

p = inputParser;
p.KeepUnmatched = true;

addRequired(p,  'phase',                         @(x) isnumeric(x));
addRequired(p,  'alpha1',                        @(x) isscalar(x));

addParameter(p, 'K',            [],              @(x) isnumeric(x));
addParameter(p, 'mu1',          defaultMu,       @(x) isscalar(x));
addParameter(p, 'mu2',          defaultMu2,      @(x) isscalar(x));
addParameter(p, 'maxOuterIter', defaultNoOuter,  @(x) isscalar(x));
addParameter(p, 'tolUpdate',    defaultTol,      @(x) isscalar(x));
addParameter(p, 'regweight',    defaultRegweight,@(x) isnumeric(x));
addParameter(p, 'weight',       defaultW,        @(x) isnumeric(x));
addParameter(p, 'magnitude',    [],              @(x) isnumeric(x));
addParameter(p, 'isPrecond',    true,            @(x) islogical(x));
addParameter(p, 'isGPU',        true,            @(x) islogical(x));

parse(p, phase, alpha, params, varargin{:});


%% -------------------------------------------------------------------------
% Parsed values
% --------------------------------------------------------------------------

phase = single(p.Results.phase);
alpha = single(p.Results.alpha1);

mu           = p.Results.mu1;
mu2          = p.Results.mu2;
maxOuterIter = p.Results.maxOuterIter;
tolUpdate    = p.Results.tolUpdate;

regweight = p.Results.regweight;

if ~any(strcmpi(p.UsingDefaults, 'regweight'))
    if ndims(regweight) == 3
        regweight = repmat(regweight, [1 1 1 3]);
    end
end

magnitude = p.Results.magnitude;

if any(strcmpi(p.UsingDefaults, 'w')) && ~any(strcmpi(p.UsingDefaults, 'magnitude'))
    W = magnitude .* magnitude;
else
    W = p.Results.weight .* p.Results.weight;
end

isGPU     = p.Results.isGPU;
isPrecond = p.Results.isPrecond;


%% -------------------------------------------------------------------------
% Dipole kernel
% --------------------------------------------------------------------------

if any(strcmpi(p.UsingDefaults, 'k'))

    if isfield(params, 'voxelSize')
        voxelSize = params.voxelSize;
    else
        voxelSize = [1 1 1];
    end

    if isfield(params, 'B0direction')
        B0direction = params.B0direction;
    else
        B0direction = [0 0 1];
    end

    Kernel = dipole_kernel_angulated(N, voxelSize, B0direction);

else

    Kernel = p.Results.K;

end

end
