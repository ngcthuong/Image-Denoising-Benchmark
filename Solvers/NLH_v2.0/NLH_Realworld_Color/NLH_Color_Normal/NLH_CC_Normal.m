function imr_final = NLH_CC_Normal(im_n)

    
   [M1,M2,M3]=size(im_n);

   
   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Colorspace in which we perform denoising. BM is applied to the first
%%%%  component and the matching information is re-used for the other two.
%%%%
if (exist('colorspace') ~= 1),
    colorspace              = 'opp'; %%% (valid colorspaces are: 'yCbCr' and 'opp')
end
   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Change colorspace, compute the l2-norms of the new color channels
%%
[zColSpace, ~] = function_rgb2LumChrom(im_n, colorspace);
im=zColSpace; 


[~, sigma_est, ~] = Noise_estimation(im_n(:,:,2)); %% estimate noisy level
sigma_est = double(sigma_est);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dis_map,~] = NL_distance(8,16,2,39,single(8),double(im_n(:,:,2)));
dis_map = double(dis_map);

%%%%%%%%%%%%%%%%%% 
Ns     = 43;
N3     = 4;
N_step = 6;
N1     = 6;
N2     = 16;
%%%%%%%%%%%%%%%%





alpha = 0.618; 

sigma_est_b = sigma_est*4.0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


    k   = 2;
    
    Thr = 2.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imp1 = im(:,:,1);
imp2 = im(:,:,2);
imp3 = im(:,:,3);
imr1 = imp1;
imr2 = imp2;
imr3 = imp3;


impp = max(im_n(:,:,2),im_n(:,:,2));

lamda = 0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% NLH_Basic_Estimation

for ll = 1:k 

        [imr1,imr2,imr3] = NLH_Color_Normal(N1,N2,N3,Ns,N_step,double(alpha*imr1+(1-alpha)*imp1),double(alpha*imr2+(1-alpha)*imp2),double(alpha*imr3+(1-alpha)*imp3),single(Thr),...
        sigma_est_b/255,double(lamda*imr1+(1-lamda)*impp));
    
        N_step = 3;
end 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% NLH_wiener

N2      = 64;
N3      = 8;
Ns      = 129;
N1      = 16;
N_step1 = 16;
N_step2 = 10;


imp1 = im(:,:,1);
imp2 = im(:,:,2);
imp3 = im(:,:,3);
 

gamma = 2.3;  
 
beta = 0.8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       [y_est1,y_est2,y_est3] = NLH_Color_Wiener_Normal(N1,N2,N3,Ns,N_step1,double(imp1),double(imp2),double(imp3),single(gamma),...
                sigma_est_b/255, double(dis_map),double(imr1*gamma),double(imr2*gamma),double(imr3*gamma),double(imr2*beta+imp1*(1-beta)));            

                        
        [y_est1,y_est2,y_est3] = NLH_Color_Wiener_Normal(N1,N2,N3,Ns,N_step2,double(y_est1),double(y_est2),double(y_est3),single(gamma),...
                sigma_est_b/255, double(dis_map),double(imr1*gamma),double(imr2*gamma),double(imr3*gamma),double(y_est1*beta+imp1*(1-beta)));
         
     
            
y_est = zeros(M1,M2,M3);
y_est(:,:,1)=y_est1(:,:);
y_est(:,:,2)=y_est2(:,:);
y_est(:,:,3)=y_est3(:,:);  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Convert back to RGB colorspace

 imr_final = function_LumChrom2rgb(y_est, colorspace);


 
 


 
 







