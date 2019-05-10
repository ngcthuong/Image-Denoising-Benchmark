Program for "Detail-preserving Image Denoising via Adaptive Clustering and Progressive PCA Thresholding".
Matlab version:R2015b or R2013a.
Author: Wenzhao Zhao

The algorithm is described in the following article:
W. Zhao, Y. Lv, Q. Liu and B. Qin, "Detail-Preserving Image Denoising via Adaptive Clustering and Progressive PCA Thresholding," in IEEE Access, vol. 6, pp. 6303-6315, 2018.
doi: 10.1109/ACCESS.2017.2780985


Contents:
demo_ACPTdenoising.m - demonstration program for image denoising
demo_NoiseLevelTest.m - demonstration program for noise level estimation

noiselevel.m - the noise variance estimation algorithm
est_patch.m
merging_thr.m
ACPT.m - ACPT algorithm
progressiveTHR.m
columns2im.m
image2cols.m
litekmeans_m.m - K-means clustering algorithm

psnr.m
ssim_index.m
FeatureSIM.m

f6.png
lena.png

-------------------------------------------------------------------
 Change log
-------------------------------------------------------------------
v1.1  (25 May 2018)
!fixed a problem with ACPT when the input is a small size image 
+Improved the usage examples and annotations

v1.0  (12 December 2017)
+Initial version