
% % % clear all
ima = double(imread('f6.png'));
level = 30;%
randn('seed',0);

rima = ima+randn(size(ima))*level;
% 
tic
    Y1 = ACPT(rima,level);
% % if noise level is unknown :
% Y1= ACPT(rima);
yt=toc;

y_PSNR = psnr(ima,Y1,255);
y_SSIM = ssim_index(ima,Y1);
y_FSIM = FeatureSIM(ima, Y1);

disp(['PSNR = ' num2str(y_PSNR) ', SSIM = ' num2str(y_SSIM) ', FSIM = ' num2str(y_FSIM)])
