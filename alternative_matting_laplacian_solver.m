%
%This package contains the code for Alternative Matting Laplacian as
% described in
%
%[Pitie16] An Alternative Matting Laplacian. F. Pitie (2016). In 
%          International Conference on Image Processing (ICCV'16)
%
% send an email to fpitie@mee.tcd.ie if you want more information
%
%
function [alpha, a, R] = alternative_matting_laplacian_solver(I, varargin)

p = inputParser;

addRequired(p,'I',@isnumeric);

addOptional(p,'trimap',      @isnumeric);
addOptional(p,'B',       [], @isnumeric);
addOptional(p,'wB',      [], @isnumeric);
addOptional(p,'F',       [], @isnumeric);
addOptional(p,'wF',      [], @isnumeric);
addOptional(p,'sigma_r',  1, @isnumeric);
addOptional(p,'T',     1e-8, @isnumeric);
addOptional(p,'alpha0',  [], @isnumeric);
addOptional(p,'walpha0', [], @isnumeric);
addOptional(p,'model',   [], @ischar);
addOptional(p,'width',   [], @isnumeric);

parse(p, I, varargin{:});

trimap    = p.Results.trimap;     
sigma_r   = p.Results.sigma_r;    
T         = p.Results.T;
walpha0   = p.Results.walpha0;
alpha0    = p.Results.alpha0;
B         = p.Results.B;
wB        = p.Results.wB;
F         = p.Results.F;
wF        = p.Results.wF;
model     = p.Results.model;
width     = p.Results.width;

if isempty(width)
    width = size(I,2);
end

resize = (width ~= size(I,2));

if isempty(model)
    model = 'affine';
end

if isempty(alpha0) && isempty(trimap)
    error('a trimap or alpha0 must be provided');
end

if isempty(alpha0)
    alpha0 = (trimap > 250) * 1 + (trimap < 5) * 0;
    walpha0 = (trimap > 250 |  trimap < 5) * 10;
end

if isempty(walpha0)
    walpha0 = 10*ones(size(alpha0));
end

fprintf('... computing local colour covariances\n');

Ihires = I;
[vres, hres, ~] = size(I);

if (strcmp(model, 'linear'))
    X = I;
elseif (strcmp(model, 'affine'))
    X = cat(3, I, ones(vres, hres));
else
    error('unknown model');
end

K    = size(X, 3);


% we denote (X X') as R
R = estimate_R(X);

% setting priors for X_B
if ~exist('wB', 'var') || isempty(wB)
    wB = 0.1 * ones(vres,hres);
end

if ~exist('B','var') || isempty(B)
    RB = zeros(size(R));
    betaB = zeros(vres, hres, K);
else
    XB = cat(3, B, ones(vres, hres));
    RB = estimate_R(XB);
    RB = wblur(RB, wB, sigma_r); 
    betaB = zeros(vres, hres, K);
end

% setting priors for X_F
if ~exist('wF', 'var') || isempty(wF)
    wF = 0.1 * ones(vres,hres);
end

if ~exist('F','var') || isempty(F)
    RF = zeros(size(R));
    betaF = zeros(vres, hres, K);
else
    XF = cat(3, F, ones(vres, hres));
    RF = estimate_R(XF);
    RF = wblur(RF, wF, sigma_r);
    betaF = wblur(XF, wF, sigma_r);
end

% setting priors on alpha
betaT(:,:,:) = X.*repmat(alpha0, [1 1 K]);
RT = wblur(R, walpha0, sigma_r);
betaT = wblur(betaT, walpha0, sigma_r);

% combining all priors
R0 =    repmat(walpha0, [1 1 K K]).* RT + ...
        repmat(wB, [1 1 K K]) .* RB + ...
        repmat(wF, [1 1 K K]) .* RF;
beta0 = repmat(walpha0, [1 1 K]).* betaT +...
        repmat(wB, [1 1 K]) .* betaB + ...
        repmat(wF, [1 1 K]) .* betaF;


if (resize)
    R = imresize(R, [NaN width ]);
    R0 = imresize(R0, [NaN width]);
    beta0 = imresize(beta0, [NaN width]);
    hres = size(beta0,2);
    vres = size(beta0,1);    
end

R = blur(R, sigma_r);

unknown = ones(vres, hres);
inds = find(unknown);

[ys, xs] = ind2sub([vres, hres], inds);
inds_{1} = sub2ind([vres, hres], ys, min(xs+1, hres));
valids{1} = xs < hres;
inds_{2} = sub2ind([vres, hres], min(ys+1, vres), xs);
valids{2} = ys < vres;
inds_{3} = sub2ind([vres, hres], ys, max(xs-1, 1));
valids{3} = xs > 1;
inds_{4} = sub2ind([vres, hres], max(ys-1, 1), xs);
valids{4} = ys > 1;

N = K*vres*hres;
Ri = reshape(R, [], K, K);
R0i = reshape(R0, [], K, K);

fprintf('... adding matting laplacian:    \n');
tic

n_all_sp_entries = 2*K*K*(sum(valids{1}) + sum(valids{2}) + sum(valids{3}) + sum(valids{4}));
all_sp_indi = zeros(n_all_sp_entries, 1);
all_sp_indj = zeros(n_all_sp_entries, 1);
all_sp_vals = zeros(n_all_sp_entries, 1);

start = 1;
for n=1:4
    % index of central point in stencil
    ind0 = inds(valids{n});
    % index of neighbouring point in stencil
    indn = inds_{n}(valids{n});
    % weight
    n_entries = sum(valids{n});
    Tij = T;
   
    for ki=1:K
        for kj=1:K
            
            Aij = (Ri(ind0, ki, kj)  + Ri(indn, ki, kj))/2 + Tij*(ki==kj);
            
            sp_ind0 = K*(ind0-1);
            sp_indn = K*(indn-1);
            
            all_sp_indi(start:start+n_entries-1) = sp_ind0+ki;
            all_sp_indj(start:start+n_entries-1) = sp_indn+kj;
            all_sp_vals(start:start+n_entries-1) = -Aij;
            start = start + n_entries;
            
            all_sp_indi(start:start+n_entries-1) = sp_ind0+ki;
            all_sp_indj(start:start+n_entries-1) = sp_ind0+kj;
            all_sp_vals(start:start+n_entries-1) = Aij;
            start = start + n_entries;
            
            fprintf('\b\b\b\b%3d%%', ceil(((ki-1)*K+kj + (n-1)*K*K)/(K*K*4)*100));
        end
    end
end

tic;
A =  sparse(all_sp_indi, all_sp_indj, all_sp_vals, N, N, n_all_sp_entries);
t = toc;

fprintf(sprintf('done in %fs \n', t));

fprintf('... adding priors:     \n');

n_all_sp_entries = hres*vres*K*K;
all_sp_indi = zeros(n_all_sp_entries, 1);
all_sp_indj = zeros(n_all_sp_entries, 1);
all_sp_vals = zeros(n_all_sp_entries, 1);

start = 1;
for ki=1:K
    for kj=1:K
        ind0 = 1:(hres*vres);
        all_sp_indi(start:start+hres*vres-1) = K*(ind0-1) + ki;
        all_sp_indj(start:start+hres*vres-1) = K*(ind0-1) + kj;
        all_sp_vals(start:start+hres*vres-1) = R0i(ind0, ki, kj);
        start = start + hres*vres;
        fprintf('\b\b\b\b%3d%%', ceil(((ki-1)*K+kj)/(K*K)*100));
    end
end
A =  A + sparse(all_sp_indi, all_sp_indj, all_sp_vals, N, N);
A = (A + A')/2;

fprintf('\n');


fprintf('... exact sparse solver\n');
tic
beta0_ = reshape(permute(beta0, [3 1 2]), [], 1);
a_ = A\beta0_;
a = permute(reshape(a_, K, vres, hres), [2 3 1]);

if (resize)
    a = imresize(a, [size(Ihires,1) size(Ihires,2)]);
end

toc
alpha = sum(a.*X, 3);

end

% Gaussian Blur
function g = blur(f, sigma)
n = ceil(3*sigma);
if (n < 1)
    g = f;
else
    x = -n:n;
    w = exp(- x.^2/2/sigma^2)/sqrt(2*pi*sigma^2);
    w = w ./ sum(w);
    g = imfilter(imfilter(f, w, 'symmetric'), w', 'symmetric');
end

end

% computes the X'X matrices
function R = estimate_R(X)
    K = size(X,3);
    R = zeros(size(X,1), size(X,2), K, K);
    eps = 1e-6;
    for i=1:K
        for j=1:K
            R(:,:,i,j) = X(:,:,i).*X(:,:,j) + (i==j)*eps;
        end
    end
end

% computes the expected value of field u given 
% confidence map w:
% namely wu = blur(u.*w)./blur(w)
function wu= wblur(u, w, sigma)

size_u = size(u);

wu = blur(u.*repmat(w, [1 1 size_u(3:end)]), sigma) ./  ...
    repmat(blur(w, sigma) + 1e-6, [1 1 size_u(3:end)]);

end

