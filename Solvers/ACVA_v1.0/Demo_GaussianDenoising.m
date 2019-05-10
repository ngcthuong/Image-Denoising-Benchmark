
clc
clear all
addpath('uitls/')
addpath('functions/')

%% The simulation experiment for Gaussian denoising
ima = double(imread('house.png'));%The range of pixel value is 0-255.

%% Add noise
level = 20;%
randn('seed',0);
rima = ima + randn(size(ima))*level;

%% Denoising
tic
    wid = 128;
    step = 32;
    RES = ACVA(rima,wid,step,level); %
%     if noise level is unknown :
%     RES = ACVA(rima,wid,step);
yt = toc;

%% Metrics
r_PSNR = psnr(ima,RES,255);
r_SSIM = ssim_index(ima,RES);
r_FSIM = FeatureSIM(ima,RES);

fprintf( 'Gaussian Denoising Results: Sigma = %4.1f, PSNR = %5.2f, SSIM = %6.4f, FSIM = %6.4f \n\n\n', level, r_PSNR, r_SSIM, r_FSIM );
