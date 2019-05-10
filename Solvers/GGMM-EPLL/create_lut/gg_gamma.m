% % Undocumented.
%
% C-A. Deledalle, S. Parameswaran, and T. Q. Nguyen, "Image
% restoration with generalized Gaussian mixture model patch
% priors", arXiv.
%
% If you use any part of this software package please cite the above
% publication
%
%
% License details as in license.txt
% ________________________________________

function ga = gg_gamma(la, nu, mode, lts)

switch mode
    case 'direct'
        ga = gg_discrepancy(0, 1, la, nu, 'mode', 'direct');
    case 'approx'
        if ~exist('lts', 'var')
            lt = load('gamma0_lookup_table.mat');
            lt.x = lt.log_la;
            lt.y = lt.log_nu;
            lt.v = lt.log_ga0;
            lt.m_x = min(lt.x);
            lt.M_x = max(lt.x);
            lt.S_x = length(lt.x);
            lt.m_y = min(lt.y);
            lt.M_y = max(lt.y);
            lt.S_y = length(lt.y);
        end
        log_ga0 = access_lt(log(la), log(nu), lt);
        ga = exp(log_ga0) + log(sqrt(2 * pi));
end
