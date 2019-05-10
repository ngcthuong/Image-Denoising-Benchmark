function h = plotpatches(x, c, varargin)
% % Function Name: plotpatches
%
%   Display a (subset of a) collection of patches
%
% Input/Output
%
%    x          a P x K real array (of gray patches)
%               the patch size P is supposed to be square.
%               if K > 400, display only a 400 subset.
%
%    c          colorline specs for the separation.
%
%    varargin   same as plotimage
%
%    h          an handle on the current axis.

% Citation:
% If you use this code please cite: 
%
% C-A. Deledalle, S. Parameswaran, and T. Q. Nguyen, "Image
% restoration with generalized Gaussian mixture model patch
% priors", arXiv.
%
% License details as in license.txt
% ________________________________________


    if ~exist('x0', 'var')
        x0 = x;
    end
    [P, K] = size(x);
    if K > 20^2
        idx = randperm(K);
        K = 20^2;
        x = x(:, idx(1:K));
    end
    N = sqrt(P);
    x = reshape(x, [N N K]);
    if sqrt(K) == floor(sqrt(K))
        b = sqrt(K);
    else
        b = ceil(sqrt(K+1));
    end
    a = ceil(K / b);
    visu = zeros(b*N, a*N);
    for k = 1:K
        visu(floor((k-1)/a)*N + (1:N), mod(k-1, a)*N + (1:N)) = ...
            x(:, :, k);
    end
    plotimage(visu, varargin{:});
    hold on
    for i = 0:b
        plot([0 a*N]+0.5, [i*N i*N] + 0.5, 'Color', c, 'LineWidth', 1);
    end
    for j = 0:a
        plot([j*N j*N] + 0.5, [0 b*N]+0.5, 'Color', c, 'LineWidth', 1);
    end
    if nargout == 0
        clear h;
    end
