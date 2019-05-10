function addpathrec(path)
% % Function Name: addpathrec
%
%   Add all directories and subdirectories recursively
%   to matlab pathdefs variable

% Citation:
% If you use this code please cite:
%
% C-A. Deledalle, S. Parameswaran, and T. Q. Nguyen, "Image
% restoration with generalized Gaussian mixture model patch
% priors", arXiv.
%
% License details as in license.txt
% ________________________________________

list = dir(path);
for k = 1:length(list)
    if ~strcmp(list(k).name(1), '.') && list(k).isdir
        addpath([path '/' list(k).name]);
        addpathrec([path '/' list(k).name]);
    end
end