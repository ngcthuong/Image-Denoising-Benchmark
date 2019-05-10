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
a2 = zeros(N, L);
b2 = zeros(N, L);

a2 = ones(N,1) * nu;
b2 = -nu .* log(la') ...
     -nu/2.*(gammaln(1./nu)-gammaln(3./nu));

close all
fancyfigure
subplot(1, 2, 1)
imagesc(log10(la), log10(nu), a2)
set(gca, 'XTickMode', 'manual');
set(gca, 'YTickMode', 'manual');
xlabel('$\lambda$');
ylabel('$\nu$');
title('$a_2$');
axis square
subplot(1, 2, 2)
imagesc(log10(la), log10(nu), b2)
set(gca, 'XTickMode', 'manual');
set(gca, 'YTickMode', 'manual');
xlabel('$\lambda$');
ylabel('$\nu$');
title('$b_2$');
axis square

return

log_la = log(la);
log_nu = log(nu);
save('a2_lookup_table_100x100.mat', 'log_la', 'log_nu', 'a2');
save('b2_lookup_table_100x100.mat', 'log_la', 'log_nu', 'b2');
