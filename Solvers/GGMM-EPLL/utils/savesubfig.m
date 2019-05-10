function savesubfig(h, dirname)
% % Function Name: savesubfig
%
%   Save figure and its subfigures to a directory
%
% Inputs:
%   h           : handle on the figure
%   direname    : name directory
%
% Description:
%   save the figure H in directory DIRNAME and all its subplots as
%     DIRNAME/main.fig
%     DIRNAME/subplot1.fig
%     DIRNAME/subplot2.fig
%     ...
%
%   for subaxes containing image data, SAVEFIGSUB will also create
%     DIRNAME/subimg1.png
%     DIRNAME/subimg2.png
%     ...

% Citation:
% If you use this code please cite:
%
% C-A. Deledalle, S. Parameswaran, and T. Q. Nguyen, "Image
% restoration with generalized Gaussian mixture model patch
% priors", arXiv.
%
% License details as in license.txt
% ________________________________________


[status, msg,~] = mkdir(dirname);
if ~status
    error(sprintf('Cannot create directory "%s": %s', dirname, msg));
end

% Prevent figure to be closed while saveing it
hcloser = get(h, 'closer');
set(h, 'closer', '');

%
savefig(h, [ dirname '/' 'main.fig' ]);
hchildren = get(h, 'children');
K = length(hchildren);
k = K;
i = 1;
while k > 0
    if strcmp(get(hchildren(k), 'type'), 'axes')
        haxe = hchildren(k);
        [img, flag] = getimage(haxe);
        switch flag
            case 0 % Not an image
            case 1 % Indexed image
                warning('case 1 not implemented yet')
            case 2 % Intensity image in standard range
                imwrite(img, [ dirname '/' sprintf('subimg%d.png', i) ]);
            case 3 % Intensity image not in standard range
                switch get(haxe, 'climmode')
                    case 'manual'
                        range = get(haxe, 'clim');
                    case 'auto'
                        range = [min(img(:)) max(img(:))];
                end
                img = (img - range(1)) / (range(2) - range(1));
                imwrite(img, [ dirname '/' sprintf('subimg%d.png', i) ]);
            case 4 % RGB Image
                imwrite(img, [ dirname '/' sprintf('subimg%d.png', i) ]);
            case 5 % Binary image
                warning('case 5 not implemented yet')
        end
        hsubfig = figure('visible', 'off');
        set(hsubfig,'CreateFcn','set(gcf,''Visible'',''on'')')
        haxe_list = haxe;
        while k - 1 > 0 && ~strcmp(get(hchildren(k - 1), 'type'), 'axes')
            switch get(hchildren(k - 1), 'type')
                case { 'legend', 'colorbar' }
                    haxe_list = [hchildren(k - 1) haxe];
                otherwise
                    warning(sprintf('cannot export "%s"', ...
                                    get(hchildren(k - 1), 'type')));
            end
            k = k - 1;
        end
        haxe_new = copyobj(haxe_list, hsubfig);
        set(haxe_new(end), 'Position', get(0, 'DefaultAxesPosition'));
        colormap(hsubfig, colormap(h));
        filename = [ dirname '/' sprintf('subplot%d.fig', i) ];
        savefig(hsubfig, filename);
        delete(hsubfig);
        i = i + 1;
    end
    k = k - 1;
end

% Restore closing property
set(h, 'closer', hcloser);
disp(sprintf('subfigures saved to "%s"', dirname));
