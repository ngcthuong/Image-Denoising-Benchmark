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

la = logspace(log10(0.001), log10(1000), L);
nu = logspace(log10(0.3), log10(2), N);

a1 = 2*ones(N, L);
b1 = estimate_b1(nu, la)';

close all
fancyfigure
subplot(1, 2, 1)
imagesc(log10(la), (nu), a1)
set(gca, 'XTickMode', 'manual');
set(gca, 'YTickMode', 'manual');
xlabel('$\log_{10} \lambda$');
ylabel('$\log_{10} \nu$');
title('$a_1$');
axis square
caxis([0.3 3]);
subplot(1, 2, 2)
imagesc(log10(la), log10(nu), b1)
set(gca, 'XTickMode', 'manual');
set(gca, 'YTickMode', 'manual');
xlabel('$\log_{10} \lambda$');
ylabel('$\log_{10} \nu$');
title('$b_1$');
axis square

return

log_la = log(la);
log_nu = log(nu);
save('a1_lookup_table_100x100.mat', 'log_la', 'log_nu', 'a1');
save('b1_lookup_table_100x100.mat', 'log_la', 'log_nu', 'b1');
