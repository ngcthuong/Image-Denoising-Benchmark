function v = access_lt(x, y, lt)
% % Function Name: access_lt
%
%
% Citation:
% If you use this code please cite:
%
% C-A. Deledalle, S. Parameswaran, and T. Q. Nguyen, "Image
% restoration with generalized Gaussian mixture model patch
% priors", arXiv.
%
% License details as in license.txt
% ________________________________________


k = (x - lt.m_x) / (lt.M_x - lt.m_x) * (lt.S_x - 1) + 1;
k(k < 1) = 1;
k(k > lt.S_x) = lt.S_x;
kf = floor(k);
kc = kf + 1;
ikf = max(kf, 1);
ikc = min(kc, lt.S_x);

l = (y - lt.m_y) / (lt.M_y - lt.m_y) * (lt.S_y - 1) + 1;
l(l < 1) = 1;
l(l > lt.S_y) = lt.S_y;
lf  = floor(l);
lc  = lf + 1;
ilf = max(lf, 1);
ilc = min(lc, lt.S_y);

v = (1 - l  + lf) .* (1 - k  + kf) .* lt.v(ikf + lt.S_x * (ilf - 1))  + ...
    (1 - l  + lf) .* (1 - kc + k)  .* lt.v(ikc + lt.S_x * (ilf - 1))  + ...
    (1 - lc + l)  .* (1 - k  + kf) .* lt.v(ikf + lt.S_x * (ilc - 1))  + ...
    (1 - lc + l)  .* (1 - kc + k)  .* lt.v(ikc + lt.S_x * (ilc - 1));
