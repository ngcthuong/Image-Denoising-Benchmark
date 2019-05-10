%--------------------------------------------------------------------------------------------------
% This is an implementation of the PGPD algorithm for image denoising.
% Author:  Jun Xu, csjunxu@comp.polyu.edu.hk
%              The Hong Kong Polytechnic University
% Please refer to the following paper if you use this code:
% Jun Xu, Lei Zhang, Wangmeng Zuo, David Zhang, and Xiangchu Feng,
% Patch Group Based Nonlocal Self-Similarity Prior Learning for Image Denoising.
% IEEE Int. Conf. Computer Vision (ICCV), Santiago, Chile, December 2015.
% Please see the file License.txt for the license governing this code.
%--------------------------------------------------------------------------------------------------
clear,clc;
nSig = 50;
% set parameters
[par, model]  =  Parameters_Setting( nSig );
% read clean image
par.I = single( imread('cameraman.png') )/255;
% generate noisy image
randn('seed',0);
par.nim =   par.I + par.nSig*randn(size(par.I));
fprintf('The initial value of PSNR = %2.4f, SSIM = %2.4f \n', csnr( par.nim*255, par.I*255, 0, 0 ),cal_ssim( par.nim*255, par.I*255, 0, 0 ));
% PGPD denoising
[im_out,par]  =  PGPD_Denoising(par,model);
% [im_out,par]  =  PGPD_Denoising_faster(par,model); % faster speed
% calculate the PSNR and SSIM
PSNR = csnr( im_out*255, par.I*255, 0, 0 );
SSIM =  cal_ssim( im_out*255, par.I*255, 0, 0 );
fprintf('Cameraman : PSNR = %2.4f, SSIM = %2.4f \n', PSNR, SSIM );