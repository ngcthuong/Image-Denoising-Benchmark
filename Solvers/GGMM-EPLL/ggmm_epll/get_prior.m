function prior_model = get_prior(name)
% % Function Name: get_prior
%
%
% Citation:
% If you use this code please cite:
%
% C-A. Deledalle, S. Parameswaran, and T. Q. Nguyen, "Image
% restoration with generalized Gaussian mixture model patch
% priors", arXiv.
%
% License details as in license.txt
% ________________________________________



switch name
    case 'gmm'
        load('data/gmm.mat', 'GS');
        for k = 1:GS.nmodels
            GS.nu{k} = 2 * ones(GS.dim, 1);
        end
        GS.nu{k} = double(GS.nu{k});
    case 'lmm'
        load('data/lmm.mat', 'GS');
    case 'hlmm_05'
        load('data/hlmm_05.mat', 'GS');
    case 'ggmm'
        load('data/ggmm.mat', 'GS');
end
prior_model.GS     = GS;
prior_model.name   = name;