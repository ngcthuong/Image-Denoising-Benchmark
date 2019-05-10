function [SigmaX, svp] = ClosedWNNM(SigmaY, C, oureps)
% Function to do the weighted WNNM 
%   + Input: 
%       - SigmaY: the diagonal eigence matrix
%       - C: Constant value 
%       - oureps: the epsilon constant 
%   + Output: 
%       - Weighted diagonal matrix : SigmaX
%       - svp:     
    temp    = (SigmaY - oureps).^2 - 4*(C - oureps*SigmaY);             % = (SigmaY + oureps)^2 - 4 * C
    ind     = find (temp>0);
    svp     = length(ind);    
    SigmaX  = max( SigmaY(ind) - oureps + sqrt(temp(ind)) , 0) / 2;
end