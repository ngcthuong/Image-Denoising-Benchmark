function h = fancylegend(varargin)
% % Function Name: fancylegend
%
%   Same as legend but with latex interpreter

% Citation:
% If you use this code please cite:
%
% C-A. Deledalle, S. Parameswaran, and T. Q. Nguyen, "Image
% restoration with generalized Gaussian mixture model patch
% priors", arXiv.
%
% License details as in license.txt
% ________________________________________


h = legend(varargin{:});
set(h, 'interpreter', 'latex');
if nargout == 0
    clear h;
end
