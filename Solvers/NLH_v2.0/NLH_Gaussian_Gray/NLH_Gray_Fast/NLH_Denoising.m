
function [y_est] = NLH_Denoising(par, mode)

    im = par.I; % im2double(imread(image_name));  %% read a noise-free image and put in intensity range [0,1]
    im1= im;
    [M,N]=size(im);
    im = par.nim;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Noise level estimation
    if par.nSig < 0
        [~, sigma_est, ~] = Noise_estimation(im); %% estimate noisy level
        sigma_est = double(sigma_est);
    else 
        sigma_est = par.nSig;
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [dis_map,~] = NL_distance(8,16,2,39,single(8),double(im));
    dis_map = double(dis_map);
    %%%%%%%%%%%%%%%%%% 
    Ns     = 43;
    N3     = 4;
    N2     = 16;
    
    %%%%%%% Main function    
    imr = im;      
    alpha = 0.618;     
    lamda = 0.8;    
    Thr = 1.45;
   
    %%%%%%%%%%%%%%%%%%% Iteration times and thresholdings
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
        if (mode == 1)
            imr = NLH_Gray_Fast(N1,N2,N3,Ns,N_step,double(alpha*imr+(1-alpha)*im),single(Thr),...
                                sigma_est, double(lamda*imr+(1-lamda)*im));
        elseif mode == 2 
            imr = NLH_Gray_Normal(N1,N2,N3,Ns,N_step,double(alpha*imr+(1-alpha)*im),single(Thr),...
                                  sigma_est, double(lamda*imr+(1-lamda)*im));
        end
         N_step = N_step -1;
         if N_step <= 3 || ll == k
             N_step = 3;
         end
    end 
    
    imr = double(imr);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% NLH_wiener
    N2 = 64;
    N3 = 8;
    Ns = 129;

    if sigma_est < 25
        N1      = 8;
        N_step1 = 8;
        N_step2 = 7;
    elseif sigma_est >= 25 && sigma_est < 75
        N1      = 16;
        N_step1 = 16;
        N_step2 = 15;    
    else
        N1      = 24;
        N_step1 = 24;
        N_step2 = 23;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta = 0.8;
    gamma = 2.3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (mode == 1)
        y_est = NLH_Gray_Fast_Wiener(N1,N2,N3,Ns,N_step1,double(im),single(gamma),...
                sigma_est/255, double(dis_map),beta,double(imr*gamma),double(imr*beta+im*(1-beta)));            
    
           
        y_est = NLH_Gray_Fast_Wiener(N1,N2,N3,Ns,N_step2,double(y_est),single(gamma),...
                sigma_est/255, double(dis_map),beta,double(imr*gamma),double(y_est*beta+im*(1-beta)));
    elseif mode == 2
         y_est = NLH_Gray_Wiener_Normal(N1,N2,N3,Ns,N_step1,double(im),single(gamma),...
                sigma_est/255, double(dis_map),beta,double(imr*gamma),double(imr*beta+im*(1-beta)));            
    
                y_est=double(y_est);
         
                    
        y_est = NLH_Gray_Wiener_Normal(N1,N2,N3,Ns,N_step2,double(y_est),single(gamma),...
                sigma_est/255, double(dis_map),beta,double(imr*gamma),double(y_est*beta+im*(1-beta)));
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









