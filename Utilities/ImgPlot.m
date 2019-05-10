% Author: Shunsuke Ono (ono@sp.ce.titech.ac.jp)

function[] = ImgPlot(x, caption, figNum, varargin)

figure(figNum)
if ~isempty(varargin)
    subplot(varargin{1}(1),varargin{1}(2),varargin{1}(3)), imshow(x, 'InitialMagnification', 'fit'), title(caption);
else
    imshow(x, 'InitialMagnification', 'fit'), title(caption);
end