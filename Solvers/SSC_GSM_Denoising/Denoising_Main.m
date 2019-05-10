% =========================================================================
% SSC-GSM-Denoising for image denoising, Version 1.0
% Copyright(c) 2015 Weisheng Dong
% All Rights Reserved.
%
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is here
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------
%
% This is an implementation of the algorithm for Hyperspectral image super-
% resolution from a pair of low-resolution hyperspectral image and a high-
% resolution RGB image.
% 
% Please cite the following paper if you use this code:
%
% W. Dong, G. Shi, Y. Ma, and X. Li, “Image Restoration via Simultaneous 
% Sparse Coding: Where Structured Sparsity Meets Gaussian Scale Mixture,” 
% International Journal of Computer Vision (IJCV), vol. 114, no. 2, 
% pp. 217-232, Sep. 2015.
% 
%--------------------------------------------------------------------------
clc;
clear;
randn('seed',0);

fn               =    'Data\Denoising_test_images\Monarch_full.tif'; 
dict             =    2; 
L                =    [5, 10, 15, 20, 50, 100];
idx              =    5;

par              =    Parameters_setting( L(idx), idx );
par.I            =    double( imread( fn ) );
par.nim          =    par.I + L(idx)*randn(size( par.I ));
    
[im PSNR SSIM]   =    SSC_GSM_Denoising( par );    

imwrite(im./255, 'Results\SSCGSM_den_Monarch.tif');
disp( sprintf('%s: PSNR = %3.2f  SSIM = %f\n', 'House', PSNR, SSIM) );



  