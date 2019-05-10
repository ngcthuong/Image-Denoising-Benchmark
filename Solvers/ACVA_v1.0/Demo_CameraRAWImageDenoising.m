
clc
clear all
addpath('uitls/')
addpath('functions/')

%% The simulation experiment for camera raw image denoising 
y = im2double(imread('lena.png'));

%% Noise parameters a,b
a = 1.0/200;   b = 0; 

%% Add noise
init = 1;
init = round(100000*rand(1));
randn('state',init);  rand('state',init);
randn('seed', init);  rand('seed', init);

% Poissonian component
chi = 1/a;
z = poissrnd(max(0,chi*y))/chi;
% Gaussian component
z = z+sqrt(b)*randn(size(y)); 

%% R, G, B --> R, G1, G2, B
[yR,yG1,yG2,yB]=sim_rawCameraRGB(y);
[zR,zG1,zG2,zB]=sim_rawCameraRGB(z);

alpha = a;
sigma = sqrt(max(b,0));
        
zz{1} = zR; 
zz{2} = zG1;
zz{3} = zG2; 
zz{4} = zB;
%% Denoising 
for zi = 1:4
    fz = GenAnscombe_forward(zz{zi}, sigma, alpha);% Apply forward generalized Anscombe variance-stabilizing transformation(http://www.cs.tut.fi/~foi/invansc/)
    
    sigma_den = 1;  % Standard-deviation value assumed after variance-stabiliation
    % Scale the image 
    scale_range = 1;
    scale_shift = (1-scale_range)/2;

    maxzans = max(fz(:));
    minzans = min(fz(:));
    fz = (fz-minzans)/(maxzans-minzans);   sigma_den = sigma_den/(maxzans-minzans);
    fz = fz*scale_range+scale_shift;       sigma_den = sigma_den*scale_range;
    
    % ACVA
    STEP_M = 128;
    STEP_m = 32;
    RES = ACVA(fz*255,STEP_M,STEP_m,sigma_den*255);  
    RES = RES/255;
    
    RES = (RES-scale_shift)/scale_range;
    RES = RES*(maxzans-minzans)+minzans;      
    
    % Inverse transformation
    z_denoi = GenAnscombe_inverse_exact_unbiased(RES,sigma,alpha);
    pro0{zi} = z_denoi;       
end

zRd = pro0{1}; 
zG1d = pro0{2};
zG2d = pro0{3};
zBd = pro0{4}; 

%% Metrics
avg_PSNR = (psnr(yR*255,zRd*255,255)+psnr(yB*255,zBd*255,255)+psnr(yG1*255,zG1d*255,255)+psnr(yG2*255,zG2d*255,255))/4.0;

avg_SSIM = (ssim_index(yR*255,zRd*255)+ssim_index(yB*255,zBd*255)+ssim_index(yG1*255,zG1d*255)+ssim_index(yG2*255,zG2d*255))/4.0;

avg_FSIM = (FeatureSIM(yR*255, zRd*255)+FeatureSIM(yB*255, zBd*255)+FeatureSIM(yG1*255, zG1d*255)+FeatureSIM(yG2*255, zG2d*255))/4.0;

fprintf( 'Simulation Results: Alpha = %6.4f, PSNR = %5.2f, SSIM = %6.4f, FSIM = %6.4f \n\n\n', alpha, avg_PSNR, avg_SSIM, avg_FSIM );

