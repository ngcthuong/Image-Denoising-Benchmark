clear;
clc;

% make sure you are in the MultiScaleEPLL directory
addpath(genpath(pwd));

% params
patchSize = 8;
noiseSD = 50/255;
betas = (1/noiseSD^2)*[1 4 8 16 32 64];

% models
load GSModel_8x8_200_2M_noDC_zeromean;
load GMM_high;
models = {GS,GS,GS};

% jump size
jmp = [1,2,4];

% filters
G = fspecial('gaussian',11,0.8);
GG = fspecial('gaussian',11,1.5);
H = zeros(size(G));
H((1+size(H,1))/2,(1+size(H,2))/2) = 1;
H = H - fspecial('gaussian',11,0.6);

filters = {1,G,GG};

% weights
weights = [1,0.05,0.05];

% load image
I = imread('Denoising_test_images\Lena512.png');
I = im2double(I);

% add noise
randn('seed',0);
noiseI = I + noiseSD*randn(size(I));

% multi-scale EPLL denoising
[cleanI,psnr] = denoise(noiseI,I,models,betas,noiseSD,jmp,filters,weights);

% present result
figure; imshow(cleanI,[]);
psnr


