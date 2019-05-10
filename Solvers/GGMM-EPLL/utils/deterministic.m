function varargout = deterministic(str, varargin)
% % Function Name: deterministic
%
%   Enter or exit deterministic code sections
%
% Inputs:
%   str         : 'on' of 'off'
%   varargin    : state to restore when str='off'
%
% Outputs:
%   varargout   : previous state when str='on'
%
% Example:
%
%   state = deterministic('on')
%
%   %%% deterministic section here%%%
%
%   deterministic('off', state)

% Citation:
% If you use this code please cite:
%
% C-A. Deledalle, S. Parameswaran, and T. Q. Nguyen, "Image
% restoration with generalized Gaussian mixture model patch
% priors", arXiv.
%
% License details as in license.txt
% ________________________________________


rng('default')
switch str
    case 'off'
        rng(varargin{1});
    otherwise
        varargout{1} = rng;
        rng(100);
end
