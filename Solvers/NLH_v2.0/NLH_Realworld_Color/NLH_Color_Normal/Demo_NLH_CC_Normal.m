
clear all
close all


image_name = [

   '5dma_iso3200_1_real.png';

   '5dma_iso3200_2_real.png';

   '5dma_iso3200_3_real.png';

   'd600_iso3200_1_real.png';

   'd600_iso3200_2_real.png';

   'd600_iso3200_3_real.png';

   'd800_iso1600_1_real.png';
 
   'd800_iso1600_2_real.png';

   'd800_iso1600_3_real.png';

   'd800_iso3200_1_real.png';

   'd800_iso3200_2_real.png';
     
   'd800_iso3200_3_real.png';

   'd800_iso6400_1_real.png';

   'd800_iso6400_2_real.png';

   'd800_iso6400_3_real.png'

    ];

im_noise  = cell(15,1);
Idenoised = cell(15,1);
im_mean   = cell(15,1);

tic,
parfor ii = 1 :15
   im_noise{ii,1} = im2double(imread(image_name(ii,:)));  %% read a noisy image
    
   Idenoised{ii,1} = NLH_CC_Normal(im_noise{ii,1});
   
end
toc,

PSNR_SUM = 0.0;
SSIM_SUM = 0.0;

for ii = 1:15
    
 im3 = im2double(imread( fullfile('/Users/houyingkun/NLH_v2.0/NLH_Realworld_Color/real_images',[image_name(ii,1:end-8),'mean.png'])));  % Changing to your own directory
 imr_final = Idenoised{ii,1};

 PSNR = psnr(im3,imr_final);
 figure,imshow(double(imr_final));
 
 PSNR_SUM = PSNR_SUM + PSNR;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate SSIM value

SSIM = ssim(im3,imr_final);

fprintf('Final Estimation, PSNR: %.2f dB, MSSIM: %.4f \n',PSNR, SSIM);  

SSIM_SUM = SSIM_SUM + SSIM;    
end

PSNR_AVE = PSNR_SUM/15
SSIM_AVE = SSIM_SUM/15

   
   
   
   
   
   
   
  
   
   
   
   
   
   
   
   
   