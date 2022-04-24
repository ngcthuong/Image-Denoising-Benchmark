
clear all
close all


image_name = [

%    'RGB_Lena512.png'
%    'RGB_House256.png'
%    'RGB_Peppers512.png'
   'RGB_barbara.bmp'
%    'RGB_Baboon512.png'
%    'RGB_F16_512.png'
%    'RGB_kodim01.png'
%    'RGB_kodim02.png'
%    'RGB_kodim03.png'
%    'RGB_kodim012.png'

    ];
    
   sigma = 50.0;

   im1 = im2double(imread(image_name));  %% read a noise-free image
    
   [M1,M2,M3]=size(im1);
   
   
    randn('seed', 0);                          %% generate seed
    im_n = im1+(sigma/255)*randn(size(im1)); %% create a noisy image
 
    figure,imshow(im_n);
   
    

if (exist('colorspace') ~= 1),
    colorspace              = 'opp'; %%% (valid colorspaces are: 'yCbCr' and 'opp')
end
   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Change colorspace, compute the l2-norms of the new color channels
%%
[zColSpace,~] = function_rgb2LumChrom(im_n, colorspace);
im=zColSpace; 


[~, sigma_est, ~] = Noise_estimation(im_n(:,:,2)); %% estimate noisy level
sigma_est = double(sigma_est);

sigma_est_b = sigma_est*1.0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dis_map,sigma_tb1] = NL_distance(8,16,2,39,single(8),double(im_n(:,:,2)));
dis_map = double(dis_map);


%%%%%%%%%%%%%%%%%% 
Ns     = 43;
N3     = 4;
N_step = 11;
N1     = 11;
N2     = 16;
%%%%%%%%%%%%%%%%
tic,

%%%%%%% Main function

alpha = 0.618;  

if sigma_est < 15
    lamda = 0.9;
elseif sigma_est >=15 && sigma_est < 20
     lamda = 0.9;
elseif sigma_est >=20 && sigma_est < 35
     lamda = 0.6;
elseif sigma_est >=35 && sigma_est < 75
     lamda = 0.6; 
else
     lamda = 0.75;
end

if sigma_est >= 75
    alpha = 0.75;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Iteration times and thresholding
if sigma_est < 7.5
    k = 2;
    Thr = 0.65;
 elseif sigma_est >= 7.5 && sigma_est < 15
    k = 3;
    Thr = 0.6; 
 elseif sigma_est >= 15 && sigma_est < 20
    k = 4;
    Thr = 0.6; 
  elseif sigma_est >= 20 && sigma_est < 30
    k = 4;
    Thr = 0.6;     
 elseif sigma_est >= 30 && sigma_est < 40
    k = 5;
    Thr = 0.6;    
 elseif sigma_est >= 40 && sigma_est < 50
    k = 5;
    Thr = 0.6;   
 elseif sigma_est >= 50 && sigma_est < 60
    k = 6;
    Thr = 0.6; 
 elseif sigma_est >= 60 && sigma_est < 75
    k = 6;
    Thr = 0.6;  
 elseif sigma_est >= 75 && sigma_est < 85
     k = 6;
     Thr = 0.4; 
 elseif sigma_est >= 85 && sigma_est < 95
     k = 7;
     Thr = 0.4;
 else  
     k = 7;
     Thr = 0.4;
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


imp1 = im(:,:,1);
imp2 = im(:,:,2);
imp3 = im(:,:,3);
imr1 = imp1;
imr2 = imp2;
imr3 = imp3;

hsv = rgb2hsv(im_n);

impp = max(imp1,im_n(:,:,2));

%%%%%%% Main function

tic,
for ll = 1:k 


        [imr1,imr2,imr3] = NLH_Gaussian_Color_Normal(N1,N2,N3,Ns,N_step,double(alpha*imr1+(1-alpha)*imp1),double(alpha*imr2+(1-alpha)*imp2),double(alpha*imr3+(1-alpha)*imp3),single(Thr),...
        sigma_est_b/255,double(lamda*imr1+(1-lamda)*impp));
       N_step = N_step -1;
        if N_step <= 3 || ll == k
           N_step = 3;
        end
end 



imr=zeros(M1,M2,M3);
imr(:,:,1)=imr1(:,:);
imr(:,:,2)=imr2(:,:);
imr(:,:,3)=imr3(:,:);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Convert back to RGB colorspace
%%%
imr_basic = function_LumChrom2rgb(imr, colorspace);



imr_basic = double(imr_basic);

PSNR = psnr(im1,imr_basic);
SSIM = ssim(im1,imr_basic);
 
fprintf( 'Basic estimation: PSNR = %2.4f, SSIM = %2.4f \n', PSNR, SSIM);
 
figure,imshow(imr_basic);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Noise level estimation

imr = rgb2gray(imr_basic);

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
elseif sigma_est >=25 && sigma_est < 75
    N1      = 16;
    N_step1 = 16;
    N_step2 = 10;   
else
    N1      = 24;
    N_step1 = 24;
    N_step2 = 20;
end



imp1 = im(:,:,1);
imp2 = im(:,:,2);
imp3 = im(:,:,3);
 

gamma = 2.0 + 10.0*sigma_est1/sigma_est;

beta = 0.8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       [y_est1,y_est2,y_est3] = NLH_Gaussian_Color_Wiener_Normal(N1,N2,N3,Ns,N_step1,double(imp1),double(imp2),double(imp3),single(gamma),...
                sigma_est_b/255, double(dis_map),double(imr1*gamma),double(imr2*gamma),double(imr3*gamma),double(imr1*beta+imp1*(1-beta)));            

                         
        [y_est1,y_est2,y_est3] = NLH_Gaussian_Color_Wiener_Normal(N1,N2,N3,Ns,N_step2,double(y_est1),double(y_est2),double(y_est3),single(gamma),...
                sigma_est_b/255, double(dis_map),double(imr1*gamma),double(imr2*gamma),double(imr3*gamma),double(y_est1*beta+imp1*(1-beta)));
         


y_est = zeros(M1,M2,M3);
y_est(:,:,1)=y_est1(:,:);
y_est(:,:,2)=y_est2(:,:);
y_est(:,:,3)=y_est3(:,:);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Convert back to RGB colorspace
%%%%
imr_final = function_LumChrom2rgb(y_est, colorspace);
figure,imshow(double(imr_final));
 
 
 PSNR = psnr(im1,imr_final);
 
 SSIM = ssim(im1,imr_final);
 
 
 fprintf( 'Final estimation: PSNR = %2.4f, SSIM = %2.4f \n', PSNR, SSIM);

 toc
 







