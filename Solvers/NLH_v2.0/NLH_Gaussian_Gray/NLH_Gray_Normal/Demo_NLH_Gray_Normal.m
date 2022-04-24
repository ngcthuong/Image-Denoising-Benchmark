
clear all
close all

% clc


image_name = [
%     'Cameraman256.png'
    'house.png'
%     'peppers256.png'
%     'montage.png'
%     'Lena512.png'
%     'barbara.png'
%     'boat.png'
%     'man.png'
%     'couple.png'
%     'hill.png' 
%     'fingerprint.png'
%     'd800_iso6400_3_gray.png'
    ];


    sigma  = 5 ; 



    im = im2double(imread(image_name));  %% read a noise-free image and put in intensity range [0,1]

    im1= im;

    [M,N]=size(im);
    
     
    randn('seed', 0);                          %% generate seed
    im = im + (sigma/255)*randn(size(im)); %% create a noisy image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 figure,imshow(im);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Noise level estimation

[~, sigma_est, ~] = Noise_estimation(im); %% estimate noisy level

sigma_est = double(sigma_est);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dis_map, ~] = NL_distance(8,16,2,39,single(8),double(im));
dis_map = double(dis_map);
%%%%%%%%%%%%%%%%%% 
Ns     = 43;
N3     = 4;
N2     = 16;

%%%%%%%%%%%%%%%%%
tic,


%%%%%%% Main function

imr = im;  

alpha = 0.618; 

lamda = 0.8;

Thr = 1.45;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Iteration times and thresholdings
  
if sigma_est < 7.5
    k = 2;
    N_step = 4;
    N1     = 9;
elseif sigma_est >= 7.5 && sigma_est < 12.5
    k = 3;
    N_step = 5;
    N1     = 9;
 elseif sigma_est >= 12.5 && sigma_est < 35
    k = 4;
    N_step = 6;
    N1     = 9;
 elseif sigma_est >= 35 && sigma_est < 55
    k = 4;
    N_step = 6;
    N1     = 9;
  elseif sigma_est >= 55 && sigma_est < 75
      k = 5;
     N_step = 7;
     N1     = 10;
 elseif sigma_est >= 75 && sigma_est < 85
     k = 6;
      N_step = 8;
      N1     = 11;
else 
      k = 7;
      N_step = 9;
      N1     = 11;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for ll = 1:k 

     
        imr = NLH_Gray_Normal(N1,N2,N3,Ns,N_step,double(alpha*imr+(1-alpha)*im),single(Thr),...
        sigma_est, double(lamda*imr+(1-lamda)*im));
     
         N_step = N_step -1;
         if N_step <= 3 || ll == k
             N_step = 3;
         end
 
end 



PSNR = psnr(im1,double(imr));  
SSIM = ssim(im1,double(imr));  

figure,imshow(imr);
 
fprintf('Basic ESTIMATE, PSNR = %2.4f, SSIM = %2.4f \n', PSNR, SSIM);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Noise level estimation

[~, sigma_est1, ~] = Noise_estimation(imr); %% estimate noisy level

sigma_est1 = double(sigma_est1);



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% NLH_wiener
N2 = 64;
N3 = 8;
Ns = 129;


if sigma_est < 25
    N1      = 8;
    N_step1 = 8;
    N_step2 = 6;
elseif sigma_est >= 25 && sigma_est < 75
    N1      = 16;
    N_step1 = 16;
    N_step2 = 10;    
else
    N1      = 20;
    N_step1 = 20;
    N_step2 = 13;
end


beta = 0.8;

gamma = 2.3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
         y_est = NLH_Gray_Wiener_Normal(N1,N2,N3,Ns,N_step1,double(im),single(gamma),...
                sigma_est/255, double(dis_map),beta,double(imr*gamma),double(imr*beta+im*(1-beta)));            
    
                y_est=double(y_est);
         
                    
        y_est = NLH_Gray_Wiener_Normal(N1,N2,N3,Ns,N_step2,double(y_est),single(gamma),...
                sigma_est/255, double(dis_map),beta,double(imr*gamma),double(y_est*beta+im*(1-beta)));
     
         
       PSNR = psnr(im1,double(y_est));        
         
         figure,imshow(y_est);
         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate SSIM value
SSIM = ssim(im1,double(y_est));

fprintf('FINAL ESTIMATE, PSNR = %2.4f, SSIM = %2.4f \n', PSNR, SSIM);



toc










