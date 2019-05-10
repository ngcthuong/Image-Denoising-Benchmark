function [ EPat, W ] = PatEstimation( NL_mat, Self_arr, Sigma_arr, CurPat, Par )
% Function to perform weighted Nuclear Norm on the set of similar patches
% + Input:
%   - NL_mat: array of index of best similar patch
%   - Self_arr: the index of current patch
%   - Sigma_arr: estimated variance of image patch
%   - CurPat: matrix of all patch
%   - Par: parameter
% + Output:
%   - Epat: estimated patch
%   - W: the weighted norm
% load test_WNNM_code
EPat = zeros(size(CurPat));
W    = zeros(size(CurPat));                                         % Weighted matrix 
for  i      =  1 : length(Self_arr)                                 % For each keypatch group
    Temp    =   CurPat(:, NL_mat(1:Par.patnum,i));                  % Non-local similar patches to the keypatch
    M_Temp  =   repmat(mean( Temp, 2 ),1,Par.patnum);               % remove mean out from image patch
    Temp    =   Temp-M_Temp;
    
    E_Temp 	=   WNNM( Temp, Par.c, Sigma_arr(Self_arr(i)), M_Temp, Par.ReWeiIter); % WNNM Estimation
    EPat(:,NL_mat(1:Par.patnum,i))  = EPat(:,NL_mat(1:Par.patnum,i)) + E_Temp;     % Add back by overlapping 
    W(:,NL_mat(1:Par.patnum,i))     = W(:,NL_mat(1:Par.patnum,i)) + ones(Par.patsize*Par.patsize,size(NL_mat(1:Par.patnum,i),1));  
end

end

