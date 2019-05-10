Program for "Texture variation adaptive image denoising with nonlocal PCA".
Matlab version:R2015b or R2013a.
Author: Wenzhao Zhao

The algorithm is described in the following article:
1) W. Zhao, Q. Liu, Y. Lv, and B. Qin, "Texture variation adaptive image denoising with nonlocal PCA," IEEE Transactions on Image Processing. 

-------------------------------------------------------------------
Contents
-------------------------------------------------------------------
Demo_GaussianDenoising.m - demonstration program for Gaussian denoising
Demo_CameraRAWImageDenoising.m - demonstration program for simulation experiments on camera raw image denoising
ACVA.m - ACVA algorithm

est_x.m
sim_rawCameraRGB.m
ProgressBar.p

psnr.m
ssim_index.m
FeatureSIM.m

lena.png
house.png

For camera raw image denoising, the following files are used to perform the generalized Anscombe variance-stabilizing transformation(http://www.cs.tut.fi/~foi/invansc/):
Anscombe_inverse_exact_unbiased.m
Anscombe_vectors.mat
GenAnscombe_forward.m
GenAnscombe_inverse_exact_unbiased.m
GenAnscombe_vectors.mat

-------------------------------------------------------------------
 Change log
-------------------------------------------------------------------
v1.0  (12 November 2018)
+Initial version