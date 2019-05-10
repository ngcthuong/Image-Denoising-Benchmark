function [ga, a1, a2, b1, b2, h, log_ga0] = gg_discrepancy_params(la, nu, lut)
% % Function Name: gg_discrepancy_params
%
%
% Look for coefficients inside the lut
% Perform bi-linear interpolation and extrapolation

% Citation:
% If you use this code please cite:
%
% C-A. Deledalle, S. Parameswaran, and T. Q. Nguyen, "Image
% restoration with generalized Gaussian mixture model patch
% priors", arXiv.
%
% License details as in license.txt
% ________________________________________


x = log(la);
y = log(nu);

linexp = true; % Enable linear extrapolation outside the lut

k = (x - lut.m_x) / (lut.M_x - lut.m_x) * (lut.S_x - 1) + 1;
if linexp
    kf = floor(k);
    kc = kf + 1;
    kf(k < 1) = 1;
    kc(k < 1) = 2;
    kf(k > lut.S_x) = lut.S_x-1;
    kc(k > lut.S_x) = lut.S_x;
else
    k(k < 1) = 1;
    k(k > lut.S_x) = lut.S_x;
    kf = floor(k);
    kc = kf + 1;
end
ikf = max(kf, 1);
ikc = min(kc, lut.S_x);

l = (y - lut.m_y) / (lut.M_y - lut.m_y) * (lut.S_y - 1) + 1;
if linexp
    lf  = floor(l);
    lc  = lf + 1;
    lf(l < 1) = 1;
    lc(l < 1) = 2;
    lf(l > lut.S_y) = lut.S_y-1;
    lc(l > lut.S_y) = lut.S_y;
else
    l(l < 1) = 1;
    l(l > lut.S_y) = lut.S_y;
    lf  = floor(l);
    lc  = lf + 1;
end
ilf = max(lf, 1);
ilc = min(lc, lut.S_y);

w1 = (1 - l  + lf) .* (1 - k  + kf);
w2 = (1 - l  + lf) .* (1 - kc + k);
w3 = (1 - lc + l)  .* (1 - k  + kf);
w4 = (1 - lc + l)  .* (1 - kc + k);
i1 = ikf + lut.S_x * (ilf - 1);
i2 = ikc + lut.S_x * (ilf - 1);
i3 = ikf + lut.S_x * (ilc - 1);
i4 = ikc + lut.S_x * (ilc - 1);

    function v = interpol(ar)
    v = w1 .* ar(i1)  + ...
        w2 .* ar(i2)  + ...
        w3 .* ar(i3)  + ...
        w4 .* ar(i4);
    end

log_ga0 = interpol(lut.log_ga0);
a1      = 2*ones(size(nu));
b1      = interpol(lut.b1);
a2      = nu;
%b2      = -nu .* log(la) - nu/2 .* (gammaln(1./nu)-gammaln(3./nu));
b2      = interpol(lut.b2);
b2(nu == 2) = b1(nu == 2);
h       = interpol(lut.h);
ga      = exp(log_ga0) + log(sqrt(2 * pi));


ga(la == 0) = log(sqrt(2 * pi));
a2(la == 0) = 2;
b1(la == 0) = -log(2);
b2(la == 0) = -log(2);
h(la == 0)  = 1;

end