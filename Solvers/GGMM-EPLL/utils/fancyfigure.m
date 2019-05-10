function h = fancyfigure(varargin)
% % Function Name: fancyfigure
%
%   Same as figure but fullsize and with latex interpreter

% Citation:
% If you use this code please cite:
%
% C-A. Deledalle, S. Parameswaran, and T. Q. Nguyen, "Image
% restoration with generalized Gaussian mixture model patch
% priors", arXiv.
%
% License details as in license.txt
% ________________________________________


h = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
set(h, 'defaulttextinterpreter', 'latex');
if nargout == 0
    clear h;
end
