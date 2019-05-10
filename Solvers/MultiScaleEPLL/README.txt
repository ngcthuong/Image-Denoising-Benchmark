--------------------------------
  Citation:
--------------------------------
This is the Matlab Package for the Multi-Scale EPLL algorithm, presented in:
Vardan Papyan and Michael Elad, Multi-Scale Patch-Based Image Restoration, IEEE Transactions on Image Processing, Vol. 25, No. 1, Pages 249-261, January 2016.

--------------------------------
  Files and their description: 
--------------------------------
demo_denoise.m - a demo of denoising
demo_deblur.m - a demo of deblurring
demo_sr.m - a demo of super-resolution

denoise.m - denoising of an image
deblur.m - deblurring of an image
sr.m - super-resolution of an image

solveInv_denoise.m - auxiliary function for denoise.m
solveInv_deblur.m - auxiliary function for deblur.m
solveInv_sr.m - auxiliary function for sr.m

my_im2col - extracts patches from an image
my_scol2im - averaging of patches extracted from an image

GSModel_8x8_200_2M_noDC_zeromean.mat - GMM model for patches extracted from natural images, taken from:
http://people.csail.mit.edu/danielzoran/

GMM_high.mat - GMM model for patches extracted from natural images filtered with a high pass filter. This was trained using the EM algorithm presented in:
http://www.mathworks.com/matlabcentral/fileexchange/26184-em-algorithm-for-gaussian-mixture-model--em-gmm-

--------------------------------
  Parameters: 
--------------------------------
Here I present the beta values used in the different applications.
The rest of the parameters can be found in the paper cited above.

################
## Denoising ##
################
1)
noiseSD = 15/255;
betas = [1 4 8 16 32];

2)
noiseSD = 25/255;
betas = [1 4 8 16 32];

3)
noiseSD = 50/255;
betas = [1 4 8 16 32 64];

4)
noiseSD = 100/255;
betas = [1 4 8 16 32 64];

################
## Deblurring ##
################
All experiments were done using the following beta values:
betas = 15*[1 2 4 8 16 32 64];

######################
## Super-resolution ##
######################
1)
psf = fspecial('gauss',11,0.6);
scale = 2;
noiseSD = 1/255;
betas = [16 32 64 128 256 512 1024 2048 4096];

2)
psf = fspecial('gauss',11,0.6);
scale = 2;
noiseSD = 5/255;
betas = [8 16 32 64 128 256 512 1024];

3)
psf = fspecial('gauss',7,1.6);
scale = 3;
noiseSD = 1/255;
betas = [4 8 16 32 64 128 256 512];

4)
psf = fspecial('gauss',7,1.6);
scale = 3;
noiseSD = 5/255;
betas = [1 2 4 8 16 32 64 128];

--------------------------------
  Contact information: 
--------------------------------
All comments are welcomed at: vardanp91@gmail.com
