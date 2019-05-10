function [E_Img]   =  WNNM_DeNoising( N_Img, O_Img, Par )
% load test.mat

E_Img           = N_Img;                                                        % Estimated Image
[Height Width]  = size(E_Img);   
TotalPatNum     = (Height-Par.patsize+1)*(Width-Par.patsize+1);                 %Total Patch Number in the image, no boundary batch
Dim             = Par.patsize*Par.patsize;  


[Neighbor_arr, Num_arr, Self_arr] =	NeighborIndex(N_Img, Par);                  % PreCompute the all the patch index in the searching window 
            NL_mat              =   zeros(Par.patnum,length(Num_arr));          % NL Patch index matrix
            CurPat              =	zeros( Dim, TotalPatNum );                  % ??
            Sigma_arr           =   zeros( 1, TotalPatNum);                     % ??
            EPat                =   zeros( size(CurPat) );                      % ??
            W                   =   zeros( size(CurPat) );                      % ??
% Org_Img patch
%[OrgPat, ~]	=	Im2Patch( O_Img, N_Img, Par );  
E_Img_All = cell(1); 
for iter = 1 : Par.Iter        
    E_Img             	=	E_Img + Par.delta*(N_Img - E_Img);                  % updated image 
    [CurPat, Sigma_arr]	=	Im2PatchWNNM( E_Img, N_Img, Par );                      % image to patch and estimate local noise variance            
    
    if (mod(iter-1,Par.Innerloop)==0)                                           % redo after a certain iteration 
        Par.patnum = Par.patnum-10;                                             % Lower Noise level, less NL patches --> reducing number NL patchs (fixed hypotehsis) --> could be improved by adaptive 
        NL_mat  =  Block_matchingWNNM(CurPat, Par, Neighbor_arr, Num_arr, Self_arr);% Caculate Non-local similar patches for each 
        if(iter==1)
            Sigma_arr = Par.nSig * ones(size(Sigma_arr));                       % First Iteration use the input noise parameter
        end
    end       

     [EPat, W]  =  PatEstimation( NL_mat, Self_arr, Sigma_arr, CurPat, Par );   % Estimate all the patches
     E_Img      =  Patch2Im( EPat, W, Par.patsize, Height, Width );             
     PSNR       = csnr( O_Img, E_Img, 0, 0 );    
    fprintf( 'Iter = %2.3f, PSNR = %2.2f \n', iter, PSNR );
       
    E_Img_All{iter} = E_Img; 
end
% file_name = [];
% save([Par.imgName '_nSig' num2str(Par.nSig) '_' num2str(Par.trial) '.mat'], 'E_Img_All'); 

return;


