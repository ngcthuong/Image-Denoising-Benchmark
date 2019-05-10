clear;
clc;

% make sure you are in the MultiScaleEPLL directory
addpath(genpath(pwd));

% params
patchSize = 8;
psf = fspecial('gauss',7,1.6);
scale = 3;
noiseSD = 5/255;
betas = [1 2 4 8 16 32 64 128];
lambda = patchSize^2/noiseSD^2;

% models
load GSModel_8x8_200_2M_noDC_zeromean;
models = {GS,GS};

% jump size
jmp = [1,2];

% filters
filters = {1,1};

% weights
weights = [1,-0.4];

% load image
dirname = 'SR_test_images/';
name = 'bike.tif';
I = double(imread([dirname name]));

% pad
pad_sz = 31;
I = padarray(I,[pad_sz,pad_sz],'symmetric');

% original image luminance channel
ori_im = rgb2ycbcr(uint8(I) );
ori_im = double(ori_im(:,:,1));

% blur and downsample
LR = Blur('fwd',I,psf);
LR = LR(1:scale:end,1:scale:end,:); 
LR = Add_noise(LR,noiseSD*255);  
B = Set_blur_matrix(scale,LR,psf);

% apply sr only on luminance channel
lr_im = rgb2ycbcr(uint8(LR));
lr_im = double(lr_im(:,:,1));

% init hr of luminance with bicubic
hr_im = imresize(lr_im,scale,'bicubic');

% init output image
[lh,lw,ch] = size(LR);
hh = lh*scale;
hw = lw*scale;
hrim = uint8(zeros(hh, hw, ch));

% bicubic
b_im = imresize(LR,scale,'bicubic');
b_im2 = rgb2ycbcr(uint8(b_im));

% copy rest channels
hrim(:,:,2) = b_im2(:,:,2);
hrim(:,:,3) = b_im2(:,:,3);

% multi-scale EPLL super-resolution
[cleanI,psnr] = sr(lr_im/255,hr_im/255,ori_im/255,models,betas,lambda,jmp,filters,weights,B,pad_sz);
psnr

% add last channel to output image
hrim(:,:,1) = uint8(cleanI*255);

% back to rgb
im_out = ycbcr2rgb(hrim);

% remove pad
im_out = im_out(1+pad_sz:end-pad_sz,1+pad_sz:end-pad_sz,:);

% present result
figure; imshow(im_out,[]);


