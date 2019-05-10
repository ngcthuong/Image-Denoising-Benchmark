function  [X] =  WNNM( Y, C, NSig, m, Iter)
% [X] =  WNNM( Y, C, NSig, m, Iter )
%   This function performs weighted nuclear norm 
%   + Input: 
%       - Y: is set of patch (removed means) 
%       - C: is a constant number = sqrt(2) 
%       - NSig: is the standard deviation the current patch 
%       - m : Mean of patch (repmat)
%       - Iter: number of iterative wieghting ??? Why not use --> set to 1
%   + Output:
%       - Estimated, thresholded patch 

    [U,SigmaY,V]  = svd(full(Y),'econ');    
    PatNum        = size(Y,2);                          % number of patch = n
    TempC         = C * sqrt(PatNum) * 2 * NSig^2;      % sqrt(2 * n) * 2 * sigma ^ 2
    [SigmaX, svp] = ClosedWNNM(SigmaY, TempC, eps);                        
    X             = U(:, 1:svp) * diag(SigmaX) * V(:,1:svp)' + m;     
return;
