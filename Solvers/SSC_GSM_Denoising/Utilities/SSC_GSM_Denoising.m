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

function [im_out, PSNR, SSIM ]   =  SSC_GSM_Denoising( par )
time0         =   clock;
nim           =   par.nim;
ori_im        =   par.I;
b             =   par.win;
[h  w ch]     =   size(nim);

N             =   h-b+1;
M             =   w-b+1;
r             =   [1:N];
c             =   [1:M]; 
disp(sprintf('PSNR of the noisy image = %f \n', csnr(nim, ori_im, 0, 0) ));

im_out      =   nim;
lamada      =   par.w;
nsig        =   par.nSig;
m           =   par.nblk;
cnt         =   1;
for iter = 1 : par.K        
    [blk_arr, wei_arr]    =   Block_matching_SSC( im_out, par);    
    
    for t = 1 : 3    
        im_out               =   im_out + lamada*(nim - im_out);
        dif                  =   im_out-nim;
        vd                   =   nsig^2-(mean(mean(dif.^2)));
        
        if (iter==1 && t==1)
            par.nSig  = sqrt(abs(vd));            
        else
            par.nSig  = sqrt(abs(vd))*par.lamada;
        end    
        X         =   Im2Patch( im_out, par );    
        Ys        =   zeros( size(X) );        
        W         =   zeros( size(X) );
        L         =   size(blk_arr,2);
        
        for  i  =  1 : L                        
            B          =   X(:, blk_arr(:, i));            
            wei        =   repmat(wei_arr(:, i)',size(B, 1), 1);
            mB         =   repmat(sum(wei.*B, 2), 1, m);
            B          =   B-mB;        
            tmp_y                    =   SSC_GSM( double(B), par.c1, par.c2, par.nSig, mB );
            Ys(:, blk_arr(1:m,i))    =   Ys(:, blk_arr(1:m,i)) + tmp_y;
            W(:, blk_arr(1:m,i))     =   W(:, blk_arr(1:m,i)) + 1;
        end

        im_out   =  zeros(h,w);
        im_wei   =  zeros(h,w);
        k        =   0;
        for i  = 1:b
            for j  = 1:b
                k    =  k+1;
                im_out(r-1+i,c-1+j)  =  im_out(r-1+i,c-1+j) + reshape( Ys(k,:)', [N M]);
                im_wei(r-1+i,c-1+j)  =  im_wei(r-1+i,c-1+j) + reshape( W(k,:)', [N M]);
            end
        end
        im_out  =  im_out./(im_wei+eps);
    
        if isfield(par,'I')
            PSNR      =  csnr( im_out, par.I, 0, 0 );
            SSIM      =  cal_ssim( im_out, par.I, 0, 0 );
%             imwrite( im_out./255, 'Results\tmp.tif' );
        end
        
        fprintf( 'Iteration %d : nSig = %2.2f, PSNR = %2.2f, SSIM = %2.4f\n', cnt, par.nSig, PSNR, SSIM );
        cnt   =  cnt + 1;
    end       
end
if isfield(par,'I')
   PSNR      =  csnr( im_out, par.I, 0, 0 );
   SSIM      =  cal_ssim( im_out, par.I, 0, 0 );
end
disp(sprintf('Total elapsed time = %f min\n', (etime(clock,time0)/60) ));
return;

function  [X W]   =   SSC_GSM( Y, c1, c2, nsig, mB )
m                  =    size(Y,2);
U                  =    getpca(Y);

Y                  =    Y(:,1:m);
A0                 =    U'*Y;
theta0             =    sqrt( sum(A0.^2, 2)/m );
theta0             =    sqrt( max( theta0.^2 - nsig^2, 0 ) );
B0                 =    (diag(1./(theta0+eps))*A0);

a                  =    sum(B0.^2, 2);
b                  =    -2*sum(B0.*A0, 2); 
c                  =    c1*nsig^2;
tmp                =    b.^2./(16*a.^2) - c./(2*a);
idx                =    tmp>=0;
tmp                =    sqrt( tmp(idx) );
a                  =    a(idx);
b                  =    b(idx);

f0                 =    c*log(eps);
t                  =    -b./(4*a);
t1                 =    t + tmp;
t2                 =    t - tmp;
f1                 =    a.*t1.^2 + b.*t1 + c*log(t1+eps);
f2                 =    a.*t2.^2 + b.*t2 + c*log(t2+eps);

ind                =    f2<f1;
f1(ind)            =    f2(ind);
t1(ind)            =    t2(ind);
ind                =    f0<f1;
t1(ind)            =    0;
theta              =    zeros( size(theta0) );
theta(idx)         =    t1;

t1                =    1./(theta.^2 + c2*nsig^2+eps);
B                 =    diag(t1)*diag(theta)'*A0;
X                 =    U*diag( theta )*B;
X                 =   (X + mB(:,1:m));
W                 =   1;
return;


