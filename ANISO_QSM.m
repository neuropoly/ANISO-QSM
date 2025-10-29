function out = ANISO_QSM(params,varargin)

totalTime = tic;

[phase, alpha, N, Kernel, mu, mu2, maxOuterIter, ...
 tolUpdate, regweight, W, ~ ,isPrecond, R2s,Drp,Drn,csf,theta] = parse_inputs_aniso(params, varargin{:});

if isfield(params,'saveMinusIter'), saveMinusIter = params.saveMinusIter; else, saveMinusIter = 25; end
saveMinusIter = max(1, min(saveMinusIter, maxOuterIter));
savedMinus = false;

% cos^2 map
c2 = cos(theta).^2;

% Weighted init for z2
Wy = (W.*phase./(W+mu2));

% Scaling (Hz/T)
B0   = 3;  gyro = 42.58;  W2 = B0*gyro;

% -------------------- Variable initialization --------------------
z1_dx = zeros(N, 'single'); z1_dy = zeros(N, 'single'); z1_dz = zeros(N, 'single');
s1_dx = zeros(N, 'single'); s1_dy = zeros(N, 'single'); s1_dz = zeros(N, 'single');

z2_dx = zeros(N, 'single'); z2_dy = zeros(N, 'single'); z2_dz = zeros(N, 'single');
s2_dx = zeros(N, 'single'); s2_dy = zeros(N, 'single'); s2_dz = zeros(N, 'single');

z3_dx = zeros(N, 'single'); z3_dy = zeros(N, 'single'); z3_dz = zeros(N, 'single');
s3_dx = zeros(N, 'single'); s3_dy = zeros(N, 'single'); s3_dz = zeros(N, 'single');

% Init from QSM 
QSM_nii = load_nii('path to QSM');
QSM = double(QSM_nii.img);
xp  = QSM;      % χ+
dX  = QSM;      % Δχ
X0  = QSM;      % χ0

% Data split variable and dual
z2 = isPrecond * Wy + (~isPrecond)*zeros(N,'single');
s2 = zeros(N,'single');

alpha_over_mu = alpha/mu;

% -------------------- Operators --------------------
[k1, k2, k3] = ndgrid(0:N(1)-1,0:N(2)-1,0:N(3)-1);
E1  = 1 - exp(2i .* pi .* k1 / N(1));
E2  = 1 - exp(2i .* pi .* k2 / N(2));
E3  = 1 - exp(2i .* pi .* k3 / N(3));
E1t = conj(E1); E2t = conj(E2); E3t = conj(E3);
EE2 = E1t .* E1 + E2t .* E2 + E3t .* E3;

fprintf('%3s\t%10s\n', 'Iter', 'Update');
tic
for t = 1:maxOuterIter
    % -------- χ+ update --------
    tx = E1t .* fftn(z1_dx - s1_dx);
    ty = E2t .* fftn(z1_dy - s1_dy);
    tz = E3t .* fftn(z1_dz - s1_dz);

    xp_prev = xp;
    xp = real(ifftn( ( ...
          mu * (tx + ty + tz) ...                                
        + mu2 * (W2) * conj(Kernel) .* fftn(z2 - s2) ...        
        + fftn( Drp .* ( Drn.*(c2.*dX + X0) + R2s ) ) ...             
    ) ./ ( eps + mu * EE2 + mu2 * (W2^2) * abs(Kernel).^2 + Drp.*Drp ) )); 
    clear tx ty tz

    % -------- Δχ update --------
    tx = E1t .* fftn(z2_dx - s2_dx);
    ty = E2t .* fftn(z2_dy - s2_dy);
    tz = E3t .* fftn(z2_dz - s2_dz);

    dX_prev = dX;
    % correct derivative w.r.t Δχ → apply c2 in spatial domain
    dX = real(ifftn( ( ...
          mu * (tx + ty + tz) ...
        + mu2 * (W2) * conj(Kernel) .* fftn( c2 .* (z2 - s2) ) ...
        + fftn( (Drn .* c2) .* ( Drp.*xp - R2s - Drn.*X0 ) ) ...
    ) ./ ( eps + mu * EE2 + mu2 * (W2^2) * abs(Kernel).^2 + Drn.*Drn )));
    clear tx ty tz

    % -------- χ0 update --------
    tx = E1t .* fftn(z3_dx - s3_dx);
    ty = E2t .* fftn(z3_dy - s3_dy);
    tz = E3t .* fftn(z3_dz - s3_dz);

    X0_prev = X0;
    X0 = real(ifftn( ( ...
          mu * (tx + ty + tz) ...
        + mu2 * (W2) * conj(Kernel) .* fftn(z2 - s2) ...
        + fftn( Drn .* ( Drp.*xp - R2s - Drn.*(c2.*dX) ) ) ...
    ) ./ ( eps + mu * EE2 + mu2 * (W2^2) * abs(Kernel).^2 + Drn.*Drn )));
    clear tx ty tz

    % ---- update ratio (not used for early break) ----
    x_update = 100 * norm([xp(:)-xp_prev(:); dX(:)-dX_prev(:); X0(:)-X0_prev(:)]) ...
                   / (norm([xp(:); dX(:); X0(:)]) + eps);
    fprintf('%3d\t%10.4f\n', t, x_update);
    if x_update < tolUpdate
        break
    end
    % -------------- TV shrinkage for χ+ -------------------------
    Fx   = fftn(xp);
    x_dx = real(ifftn(E1 .* Fx));
    x_dy = real(ifftn(E2 .* Fx));
    x_dz = real(ifftn(E3 .* Fx));
    z1_dx = max(abs(x_dx + s1_dx) - regweight(:,:,:,1)*alpha_over_mu, 0) .* sign(x_dx + s1_dx);
    z1_dy = max(abs(x_dy + s1_dy) - regweight(:,:,:,2)*alpha_over_mu, 0) .* sign(x_dy + s1_dy);
    z1_dz = max(abs(x_dz + s1_dz) - regweight(:,:,:,3)*alpha_over_mu, 0) .* sign(x_dz + s1_dz);
    s1_dx = s1_dx + x_dx - z1_dx;  s1_dy = s1_dy + x_dy - z1_dy;  s1_dz = s1_dz + x_dz - z1_dz;
    clear x_dx x_dy x_dz

    % -------------- TV shrinkage for Δχ -------------------------
    Fx   = fftn(dX);
    x_dx = real(ifftn(E1 .* Fx));
    x_dy = real(ifftn(E2 .* Fx));
    x_dz = real(ifftn(E3 .* Fx));
    z2_dx = max(abs(x_dx + s2_dx) - regweight(:,:,:,1)*alpha_over_mu, 0) .* sign(x_dx + s2_dx);
    z2_dy = max(abs(x_dy + s2_dy) - regweight(:,:,:,2)*alpha_over_mu, 0) .* sign(x_dy + s2_dy);
    z2_dz = max(abs(x_dz + s2_dz) - regweight(:,:,:,3)*alpha_over_mu, 0) .* sign(x_dz + s2_dz);
    s2_dx = s2_dx + x_dx - z2_dx;  s2_dy = s2_dy + x_dy - z2_dy;  s2_dz = s2_dz + x_dz - z2_dz;
    clear x_dx x_dy x_dz

    % -------------- TV shrinkage for χ0 -------------------------
    Fx   = fftn(X0);
    x_dx = real(ifftn(E1 .* Fx));
    x_dy = real(ifftn(E2 .* Fx));
    x_dz = real(ifftn(E3 .* Fx));
    z3_dx = max(abs(x_dx + s3_dx) - regweight(:,:,:,1)*alpha_over_mu, 0) .* sign(x_dx + s3_dx);
    z3_dy = max(abs(x_dy + s3_dy) - regweight(:,:,:,2)*alpha_over_mu, 0) .* sign(x_dy + s3_dy);
    z3_dz = max(abs(x_dz + s3_dz) - regweight(:,:,:,3)*alpha_over_mu, 0) .* sign(x_dz + s3_dz);
    s3_dx = s3_dx + x_dx - z3_dx;  s3_dy = s3_dy + x_dy - z3_dy;  s3_dz = s3_dz + x_dz - z3_dz;
    clear x_dx x_dy x_dz Fx

    % -------------- z2 data-consistency (nonlinear) --------------
    chi_minus = c2.*dX + X0;
    pred_total = real(ifftn(Kernel .* fftn(xp + chi_minus)));
    numerator   = W .* sin(z2 - phase) + mu2 .* z2 - mu2 .* (pred_total + s2);
    denominator = W .* cos(z2 - phase) + mu2 + eps;
    z2 = z2 - numerator ./ denominator;

    s2 = s2 + pred_total - z2;

    % -------------- update Drp, Drn ------------------------------
    A = [xp(:), -(chi_minus(:))]; 
    coef = (A'*A + 1e-8*eye(2)) \ (A'*R2s(:));
    Drp = max(coef(1), 0);
    Drn = max(coef(2), 0);

    
if (~savedMinus) && (t >= saveMinusIter)
    csf_mask = csf > 0;
    chi_minus = c2.*dX + X0;          
    sum_map   = xp + chi_minus;

    if any(csf_mask(:))
        beta_minus = mean(sum_map(csf_mask), 'all');
    else
        beta_minus = mean(sum_map(:), 'all');
    end

    chi_minus_early = chi_minus - beta_minus/2;
    chi_minus_early(chi_minus_early > 0) = 0;

    out.x_minus_early     = chi_minus_early;  
    out.beta_minus        = beta_minus;
    out.iter_minus_saved  = t;
    savedMinus = true;
    
end

end

% -------------------- Post-processing / outputs --------------------
chi_minus = c2.*dX + X0;
csf_mask  = csf > 0;
sum_map   = xp + chi_minus;

if any(csf_mask(:))
    beta = mean(sum_map(csf_mask), 'all');
else
    beta = mean(sum_map(:), 'all');
end

xp = xp - beta/2;
chi_minus = chi_minus - beta/2;

xp(xp<0) = 0;
chi_minus(chi_minus>0) = 0;

out.time      = toc; toc
out.totalTime = toc(totalTime);
out.x_plus    = xp;           
out.deltaX    = dX;
out.X0        = X0;
out.x_minus   = chi_minus;    
out.iter      = t;
out.beta_final = beta;

end

% ===================== INPUT PARSER (ANISOTROPY) =====================
function [phase, alpha, N, Kernel, mu, mu2, maxOuterIter, tolUpdate, regweight, W, isGPU, isPrecond, R2s,Drp,Drn,csf,theta] = parse_inputs_aniso(params, varargin)
if isstruct(params)
    if isfield(params,'input'),   phase = params.input; else, error('Provide params.input'); end
    if isfield(params,'alpha1'),  alpha = params.alpha1; else, error('Provide params.alpha1'); end
    if isfield(params,'R2s'),     R2s   = params.R2s;   else, error('Provide R2s map.'); end
    if isfield(params,'csf'),     csf   = params.csf;   else, error('Provide csf mask.'); end
    if isfield(params,'Drp'),     Drp   = params.Drp;   else, error('Provide scalar Drp.'); end
    if isfield(params,'Drn'),     Drn   = params.Drn;   else, error('Provide scalar Drn.'); end
    if isfield(params,'theta'),   theta = params.theta; else, error('Provide theta map (radians).'); end
else
    error('Provide a struct params with required fields.');
end
N = size(phase);

defaultMu  = 100*alpha;
defaultMu2 = 1.0;
defaultNoOuter = 150;
defaultTol = 0.1;
defaultRegweight = ones([N 3]);
defaultW = ones(N);

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

phase = single(p.Results.phase);
alpha = single(p.Results.alpha1);

if any(strcmpi(p.UsingDefaults, 'k'))
    if isfield(params,'voxelSize'),   voxelSize = params.voxelSize; else, voxelSize = [1,1,1]; end
    if isfield(params,'B0direction'), B0direction = params.B0direction; else, B0direction = [0,0,1]; end
    Kernel = dipole_kernel_angulated( N, voxelSize, B0direction );
else
    Kernel = p.Results.K;
end

mu           = p.Results.mu1;
mu2          = p.Results.mu2;
maxOuterIter = p.Results.maxOuterIter;
tolUpdate    = p.Results.tolUpdate;

regweight = p.Results.regweight;
if ~any(strcmpi(p.UsingDefaults, 'regweight'))
    if ndims(regweight) == 3
        regweight = repmat(regweight,[1,1,1,3]);
    end
end

magnitude = p.Results.magnitude;
if any(strcmpi(p.UsingDefaults, 'w')) && ~any(strcmpi(p.UsingDefaults, 'magnitude'))
    W = magnitude.*magnitude;
else
    W = p.Results.weight.*p.Results.weight;
end

isGPU     = p.Results.isGPU;
isPrecond = p.Results.isPrecond;
end
