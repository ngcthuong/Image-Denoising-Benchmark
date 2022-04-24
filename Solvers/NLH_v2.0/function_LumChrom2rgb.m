function yRGB=function_LumChrom2rgb(x,colormode)
% Inverse color-space transformation   ( forward transformation is function_rgb2LumChrom.m )
%
% Alessandro Foi - Tampere University of Technology - 2005 - 2006   Public release v1.03 (March 2006)
% -----------------------------------------------------------------------------------------------------------------------------------------------
%
% SYNTAX:
%
%   yRGB = function_LumChrom2rgb(x,colormode);
%
% INPUTS:
%  x  is color-transformed image (with range typically included in or equal to [0 1]^3, depending on the transformation matrix)
%
%  colormode = 'opp', 'yCbCr', or a custom 3x3 matrix (e.g. provided by the forward transform when 'pca' is selected)
%
%       'opp'      opponent color space ('opp' is equirange version)
%       'yCbCr'    standard yCbCr (e.g. for JPEG images)
%
% OUTPUTS:
%   x  is RGB image (with range [0 1]^3)
%
%
% NOTE:    'opp' is used by default if no colormode is specified
%

if nargin==1
    colormode='opp';
end
if size(colormode)==[3 3]
    A=colormode;
    B=inv(A);
else
    if strcmp(colormode,'opp')
        A =[1/3 1/3 1/3; 0.5  0  -0.5; 0.25  -0.5  0.25];
        B =[1 1 2/3;1 0 -4/3;1 -1 2/3];
    end
    if strcmp(colormode,'yCbCr')
        A=[0.299   0.587   0.114;   -0.16873660714285  -0.33126339285715   0.5;   0.5  -0.4186875  -0.0813125];
        B=inv(A);
    end
end

%%%% Make sure that each channel's intensity range is [0,1]
maxV = sum(A.*(A>0),2);
minV = sum(A.*(A<0),2);
xNormal = reshape(x,[size(x,1)*size(x,2) 3]) * diag(maxV-minV) +  repmat(minV, [1 size(x,1)*size(x,2)])'; % put in range [0,1]
yRGB = reshape(xNormal * B', [ size(x,1) size(x,2) 3]);

return;