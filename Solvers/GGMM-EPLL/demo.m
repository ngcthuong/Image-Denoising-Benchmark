% % Denoising demo script as explained in:
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
addpathrec('.')
deterministic('on');

% Parameters
sig    = 30;

% Load and generate images
x      = double(imread('cameraman.png'))/255;
[M, N] = size(x);
sig    = sig / 255;
y      = x + sig * randn(M, N);

% Load prior computed offline
prior_model{1} = get_prior('gmm');
prior_model{2} = get_prior('lmm');
prior_model{3} = get_prior('ggmm');

% Run GGMM EPLL
for k = 1:length(prior_model)
    tstart = tic;
    xhat{k} = ggmm_epll(y, sig, prior_model{k});
    toc(tstart);
end

% Display
fancyfigure;
subplot(2,2,1)
plotimage(y, 'range', [0 1]);
title(sprintf('Noisy image (PSNR %.2f, SSIM %.3f)', ...
              psnr(y, x), ...
              ssim(y, x)));
for k = 1:length(prior_model)
    subplot(2,2,1+k)
    plotimage(xhat{k}, 'range', [0 1]);
    title(sprintf('%s+EPLL (PSNR %.2f, SSIM %.3f)', ...
                  upper(prior_model{k}.name), ...
                  psnr(xhat{k}, x), ...
                  ssim(xhat{k}, x)));
end
linkaxes
