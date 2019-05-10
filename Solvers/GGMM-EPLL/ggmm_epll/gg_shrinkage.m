function x = gg_shrinkage(y, si, la, nu, mode)
% % Function Name: gg_discrepancy
%
%
% Inputs:
%   y           : real value or array
%   si          : Gaussian noise std
%   la          : GGD standard deviation
%   nu          : Shape parameter of the GGD
%   mode='trunc': trunc nu to neirest neightbor [.3 1 4/3 3/2 2 3]
%
% Outputs:
%   f           : value of the discrepancy function

% Citation:
% If you use this code please cite:
%
% C-A. Deledalle, S. Parameswaran, and T. Q. Nguyen, "Image
% restoration with generalized Gaussian mixture model patch
% priors", arXiv.
%
% License details as in license.txt
% ________________________________________


% Retrieve arguments
if ~exist('mode', 'var')
    mode = 'trunc';
end
[K, n] = size(y);
if isscalar(nu)
    nu = nu * ones(K, 1);
end
if isscalar(la)
    la = la * ones(K, 1);
end
if isscalar(si)
    si = si * ones(K, 1);
end

%% Approximation for speed up
if strcmp(mode, 'trunc')
    nu(nu < 0.3)                     = 0.3;
    nu(1        < nu & nu <= 1+2/12) = 1;
    nu(1 + 2/12 < nu & nu <= 1+5/12) = 1.33;
    nu(1 + 5/12 < nu & nu <= 1+9/12) = 1.5;
    nu(1 + 9/12 < nu & nu < 2.5)     = 2;
    nu(2.5 < nu)                     = 3;
end

% Define l2 and l1 shrinkage
l2_shrinkage = @(y,si,la) ...
    bsxfun(@times, la.^2 ./ (la.^2 + si.^2),  y);
l1_shrinkage = @(y,si,la) ...
    sign(y) .* max(bsxfun(@minus, abs(y), sqrt(2)*si.^2./la), 0);

% Core
x = zeros(K, n);

mask = true * ones(size(nu));

idx = nu < 1 & mask;
if sum(idx) > 0
    x(idx, :) = llt1_shrinkage(y(idx, :), si(idx), la(idx), nu(idx));
end
mask(idx) = false;

idx = nu == 1 & mask;
if sum(idx) > 0
    x(idx, :) = l1_shrinkage(y(idx, :), si(idx), la(idx));
end
mask(idx) = false;

idx = nu == 1.33 & mask;
if sum(idx) > 0
    x(idx, :) = l1p33_shrinkage(y(idx, :), si(idx), la(idx), nu(idx));
end
mask(idx) = false;

idx = nu == 1.5 & mask;
if sum(idx) > 0
    x(idx, :) = l1p5_shrinkage(y(idx, :), si(idx), la(idx), nu(idx));
end
mask(idx) = false;

idx = nu == 2 & mask;
if sum(idx) > 0
    x(idx, :) = l2_shrinkage(y(idx, :), si(idx), la(idx));
end
mask(idx) = false;

idx = nu == 3 & mask;
if sum(idx) > 0
    x(idx, :) = l3_shrinkage(y(idx, :), si(idx), la(idx), nu(idx));
end
mask(idx) = false;

idx = 1 < nu & nu < 2 & mask;
if sum(idx) > 0
    x(idx, :) = lin12_shrinkage(y(idx, :), si(idx), la(idx), nu(idx));
end
mask(idx) = false;

idx = 2 < nu & mask;
if sum(idx) > 0
    x(idx, :) = lgt2_shrinkage(y(idx, :), si(idx), la(idx), nu(idx));
end
mask(idx) = false;

x = real(x);

% l_nu norm with nu < 1
function x = llt1_shrinkage(y, si, la, nu)

etaln = 1/2 * (gammaln(3./nu) - gammaln(1./nu));
C = log(2-nu) - (1-nu)./(2-nu) .* log(2-2*nu) + nu./(2-nu) .* etaln;
a = nu .* exp(nu .* etaln);
t = exp(C) .* si.^(2./(2-nu)) .* la.^(-nu./(2-nu));
x = sign(y) .* ...
    bsxfun(@power, abs(y), nu-1) .* ...
    bsxfun(@minus, ...
           bsxfun(@power, abs(y), 2-nu), ...
           a .* si.^2 .* la.^(-nu));
x(bsxfun(@lt, abs(y), t)) = 0;

% l_nu norm with nu = 1.33...
function x = l1p33_shrinkage(y, si, la, nu)

etaln = 1/2 * (gammaln(3./nu) - gammaln(1./nu));
a = exp(2 * log(si) - nu .* log(la) + log(nu) + nu .* etaln);
e = sqrt(y.^2  + 4 * a.^3 / 27);
x = y + ...
    a ./ 2^(1/3) .* ...
    ((e - y).^(1/3) - (e + y).^(1/3));

% l_nu norm with nu = 1.5
function x = l1p5_shrinkage(y, si, la, nu)

etaln = 1/2 * (gammaln(3./nu) - gammaln(1./nu));
a = exp(2 * log(si) - nu .* log(la) + log(nu) + nu .* etaln);
x = sign(y) .* ...
    1/4 .* (bsxfun(@plus, -a, ...
                 sqrt(bsxfun(@plus, a.^2, 4 * abs(y))))).^2;

% l_nu norm with 1 < nu < 2
function x = lin12_shrinkage(y, si, la, nu)

[K, n]    = size(y);
etaln     = 1/2 * (gammaln(3./nu) - gammaln(1./nu));
alpha     = exp(2 * log(si) - nu .* log(la) + log(nu) + nu .* etaln);
nu        = nu * ones(1, n);

l1thres   = sqrt(2) * si.^2 ./ la;
idx       = bsxfun(@lt, abs(y), l1thres);
q         = ones(K, n);
q(idx)    = 1 ./ (nu(idx) - 1);

x         = sign(y) .* halley(q, q .* (nu-1), alpha, abs(y)).^(q);

% l_nu norm with 2 < nu
function x = lgt2_shrinkage(y, si, la, nu)

[K, n]    = size(y);
etaln     = 1/2 * (gammaln(3./nu) - gammaln(1./nu));
alpha     = exp(2 * log(si) - nu .* log(la) + log(nu) + nu .* etaln);
nu        = nu * ones(1, n);

l1thres   = sqrt(2) * si.^2 ./ la;
idx       = bsxfun(@gt, abs(y), l1thres);
q         = ones(K, n);
q(idx)    = 1 ./ (nu(idx) - 1);

x         = sign(y) .* halley(q, q .* (nu-1), alpha, abs(y)).^(q);


% l_nu norm with nu = 3
function x = l3_shrinkage(y, si, la, nu)

etaln = 1/2 * (gammaln(3./nu) - gammaln(1./nu));
a = exp(2 * log(si) - nu .* log(la) + log(nu) + nu .* etaln);
x = sign(y) .* ...
   bsxfun(@rdivide, -1 + sqrt(1 + 4 * bsxfun(@times, a, abs(y))), 2 * a);


% Halley's root finding methods
function z = halley(q, p, a, tau)

f   = @(z) z.^q + bsxfun(@times, a, z.^p) - tau;
fp  = @(z) q .* z.^(q-1) + bsxfun(@times, a, p .* z.^(p-1));
fpp = @(z) q .* (q-1) .* z.^(q-2) + bsxfun(@times, a, p .* (p - 1) .* z.^(p-2));

z = tau;
for k = 1:10
    z = z - 2 * f(z) .* fp(z) ./ (2 * fp(z).^2 - f(z) .* fpp(z));
end
idx = isnan(z);
z(idx) = 0;
