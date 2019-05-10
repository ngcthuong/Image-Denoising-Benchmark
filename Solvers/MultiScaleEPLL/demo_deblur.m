clear;
clc;

% make sure you are in the MultiScaleEPLL directory
addpath(genpath(pwd));

% params
patchSize = 8;
% K = ones(9);
% K = K./sum(K(:));
K = fspecial('gaussian',21,1.6);
noiseSD = sqrt(2)/255;
betas = 15*[1 2 4 8 16 32 64];

% models
load GSModel_8x8_200_2M_noDC_zeromean;
models = {GS,GS,GS};

% jump size
jmp = [1,2,4];

% filters
filters = {1,1,1};

% weights
weights = [1,-0.25,0.02];

% load image
dirname = 'Deblurring_test_images/';
name = 'barbara.tif';
I = imread([dirname name]);
if size(I,3)>1
    I = rgb2gray(I);
end

% pad
I = padarray(I,size(K),'symmetric');

I = double(I)/255;

% convolve with kernel and add noise
ks = floor((size(K, 1) - 1)/2);
y = conv2(I, K, 'valid');
randn('seed',0);
y = y + noiseSD*randn(size(y));
y = double(uint8(y .* 255))./255;

% code excerpt taken from Krishnan et al.

% edgetaper to better handle circular boundary conditions
y = padarray(y, [1 1]*ks, 'replicate', 'both');
for a=1:4
  y = edgetaper(y, K);
end

noiseI = y;

% multi-scale EPLL deblurring
[cleanI,psnr] = deblur(noiseI,I,K,models,betas,noiseSD,jmp,filters,weights);

% present result
figure; imshow(cleanI(1+size(K,1):end-size(K,1),1+size(K,2):end-size(K,2)),[]);
psnr


