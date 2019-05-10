function h = plotimage(img, varargin)
% % Function Name: plotimage
%
%   Display an image
%
% Inputs:
%   img         : a M x N array
%
% Outputs:
%   h           : handle on the created axes
%
% Optional arguments:
%   range       : range on wich pixel values are mapped to grey or
%                 RGB values (default [0 255])
%
%   alpha       : for gamma correction: beta * img.^alpha (default 1)
%
%   beta        : for gamma correction: beta * img.^alpha (default 1)
%
%   Q           : number of auqntification levels (default 256).
%                 incompatible with on Windows systems
%
%   adjust      : redefine the range automatically:
%                 'no':   keep the original range
%                 'auto': range is determined form extreme values
%                 'm1s':  range is [0 m+s] where m and s are the mean
%                         and standard deviation of img
%                 'm2s':  range is [0 m+2*s]
%                 'm3s':  range is [0 m+3*s]
%                 'm4s':  range is [0 m+4*s]
%                 'm5s':  range is [0 m+5*s]
%
%   rangeof     : handle to an existing axes (created by plotimage*)
%                 used to display img with the exact same color palette.
%                 all other optional arguments are ignored.

% Citation:
% If you use this code please cite:
%
% C-A. Deledalle, S. Parameswaran, and T. Q. Nguyen, "Image
% restoration with generalized Gaussian mixture model patch
% priors", arXiv.
%
% License details as in license.txt
% ________________________________________



  if nargin < 2
    options = struct([]);
  end

  import tools.*

  options   = makeoptions(varargin{:});
  adjust    = 'no';
  if isfield(options, 'rangeof')
      h     = getoptions(options, 'rangeof', [], true);
      range = getfield(get(h, 'UserData'), 'range');
      alpha = getfield(get(h, 'UserData'), 'alpha');
      beta  = getfield(get(h, 'UserData'), 'beta');
      Q     = getfield(get(h, 'UserData'), 'Q');
  else
      range = getoptions(options, 'range', [0 255]);
      if ~isfield(options, 'range')
          adjust = getoptions(options, 'Adjust', 'no');
      end
      alpha = getoptions(options, 'alpha', 1);
      beta  = getoptions(options, 'beta', 1);
      Q     = getoptions(options, 'Q', 256);
  end
  m = range(1);
  M = range(2);

  img = beta * img.^alpha;
  switch adjust
      case 'no'
      case 'auto'
          m = min(img(:));
          M = max(img(:));
      case 'm5s'
          m = 0;
          M = mean(img(:)) + 5*std(img(:));
      case 'm4s'
          m = 0;
          M = mean(img(:)) + 4*std(img(:));
      case 'm3s'
          m = 0;
          M = mean(img(:)) + 3*std(img(:));
      case 'm2s'
          m = 0;
          M = mean(img(:)) + 2*std(img(:));
      case 'm1s'
          m = 0;
          M = mean(img(:)) + 1*std(img(:));
      case 'm0s'
          m = 0;
          M = mean(img(:)) + 0*std(img(:));
      otherwise
          error(sprintf('Display adjustment %s unknown', adjust));
  end
  if m >= M
      m = -1;
      M = 1;
  end
  if size(img, 3) == 1
      h = imagesc(img);
      caxis([m M]);
      colormap(gray(Q));
  else
      img(img > M) = M;
      img(img < m) = m;
      h = imshow((img - m) / (M - m), [0 1]);
  end
  set(h, 'UserData', struct('range', [m, M], ...
                            'alpha', alpha, ...
                            'beta', beta, ...
                            'Q', Q));
  axis image;
  axis off;

  if nargout == 0
      clear h;
  end
