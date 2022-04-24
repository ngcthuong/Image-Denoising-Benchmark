function [y, A, l2normLumChrom]=function_rgb2LumChrom(xRGB, colormode)
% Forward color-space transformation   ( inverse transformation is function_LumChrom2rgb.m )
%
% Alessandro Foi - Tampere University of Technology - 2005 - 2006   Public release v1.03 (March 2006)
% -----------------------------------------------------------------------------------------------------------------------------------------------
%
% SYNTAX:
%
%   [y A l2normLumChrom] = function_rgb2LumChrom(xRGB, colormode);
%
% INPUTS:
%   xRGB  is RGB image with range [0 1]^3
%
%   colormode = 'opp', 'yCbCr', 'pca', or a custom 3x3 matrix
%
%       'opp'     Opponent color space ('opp' is equirange version)
%       'yCbCr'   The standard yCbCr (e.g. for JPEG images)
%       'pca'     Principal components   (note that this transformation is renormalized to be equirange) 
%
% OUTPUTS:
%   y  is color-transformed image (with range typically included in or equal to [0 1]^3, depending on the transformation matrix)
%
%   l2normLumChrom (optional) l2-norm of the transformation (useful for noise std calculation)
%   A  transformation matrix  (used necessarily if colormode='pca')
%
%   NOTES:  -  If only two outputs are used, then the second output is l2normLumChrom, unless colormode='pca';
%           -  'opp' is used by default if no colormode is specified.
%
%
% USAGE EXAMPLE FOR PCA TRANSFORMATION:
%  %%%%  -- forward color transformation --
%    if colormode=='pca'
%       [zLumChrom colormode] = function_rgb2LumChrom(zRGB,colormode); % 'colormode' is assigned a 3x3 transform matrix
%    else
%       zLumChrom = function_rgb2LumChrom(zRGB,colormode);
%    end
%
%  %%%% [ ... ]  Some processing  [ ... ]
%
%  %%%%  -- inverse color transformation --
%    zRGB=function_LumChrom2rgb(zLumChrom,colormode);
%

if nargin==1
    colormode='opp';
end
change_output=0;
if size(colormode)==[3 3]
    A=colormode;
    l2normLumChrom=sqrt(sum(A.^2,2));
else
    if strcmp(colormode,'opp')
        A=[1/3 1/3 1/3; 0.5  0  -0.5; 0.25  -0.5  0.25];
    end
    if strcmp(colormode,'yCbCr')
        A=[0.299   0.587   0.114;   -0.16873660714285  -0.33126339285715   0.5;   0.5  -0.4186875  -0.0813125];
    end
    if strcmp(colormode,'pca')
        A=princomp(reshape(xRGB,[size(xRGB,1)*size(xRGB,2) 3]))';
        A=A./repmat(sum(A.*(A>0),2)-sum(A.*(A<0),2),[1 3]);  %% ranges are normalized to unitary length;
    else
        if nargout==2
            change_output=1;
        end
    end
end

%%%% Make sure that each channel's intensity range is [0,1]
maxV = sum(A.*(A>0),2);
minV = sum(A.*(A<0),2);
yNormal = (reshape(xRGB,[size(xRGB,1)*size(xRGB,2) 3]) * A' - repmat(minV, [1 size(xRGB,1)*size(xRGB,2)])') * diag(1./(maxV-minV)); % put in range [0,1]
y = reshape(yNormal, [size(xRGB,1) size(xRGB,2) 3]);

%%%% The l2-norm of each of the 3 transform basis elements 
l2normLumChrom = diag(1./(maxV-minV))*sqrt(sum(A.^2,2));

if change_output
    A=l2normLumChrom;
end

return;