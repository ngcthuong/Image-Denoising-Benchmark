function options = makeoptions(varargin)
% % Function Name: makeoptions
%
%   Extracts optional arguments
%
% Input:
%   varargin    : an even number of arguments
%                 makeoptions('arg1', value1, 'arg2', value2, ...)
%
% Output:
%   options     : a structure such that
%                 options.arg1 = value1
%                 options.arg2 = value2
%                 ...

% Citation:
% If you use this code please cite:
%
% C-A. Deledalle, S. Parameswaran, and T. Q. Nguyen, "Image
% restoration with generalized Gaussian mixture model patch
% priors", arXiv.
%
% License details as in license.txt
% ________________________________________


if mod(length(varargin), 2) == 1
    error('the number of options should be even');
end
options = struct();
rev = @(x) x(end:-1:1);
for k = rev(1:2:length(varargin))
    options = setfield(options, lower(varargin{floor(k/2)*2+1}), ...
                                varargin{floor(k/2)*2+2});
end
