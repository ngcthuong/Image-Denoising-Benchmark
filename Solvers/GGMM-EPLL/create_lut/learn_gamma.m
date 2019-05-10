% % Undocumented.
%
% C-A. Deledalle, S. Parameswaran, and T. Q. Nguyen, "Image
% restoration with generalized Gaussian mixture model patch
% priors", arXiv.
%
% If you use any part of this software package please cite the above
% publication
%
%
% License details as in license.txt
% ________________________________________

clear all
close all

cd(fileparts(mfilename('fullpath')));
addpath('..');
addpathrec('..');

L = 100;
N = 100;

la = [logspace(log10(0.001), log10(1000), L)];
nu = logspace(log10(0.3), log10(2), N);
ga = zeros(N, L);
for k = 1:length(nu)
    ga(k,:) = gg_discrepancy(0, 1, la, nu(k), 'mode', 'direct') - log(sqrt(2 * pi));
end

fancyfigure
surf(log10(la), log10(nu), log(ga));
xlabel('$\lambda$');
ylabel('$\nu$');
zlabel('$\log \gamma(\lambda, \nu)$');
xlim([-3 3]);
ylim([log10(0.3) log10(4)]);

return
log_la = log(la);
log_nu = log(nu);
log_ga0 = log(ga)';
save('gamma0_lookup_table_100x100.mat', 'log_la', 'log_nu', 'log_ga0');
