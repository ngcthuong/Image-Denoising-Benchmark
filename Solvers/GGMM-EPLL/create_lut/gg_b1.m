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

function b1 = gg_b1(la, nu, lts)

    if ~exist('lts', 'var')
        lt = load('b1_lookup_table.mat');
        lt.x = lt.log_la;
        lt.y = lt.log_nu;
        lt.v = lt.b1;
        lt.m_x = min(lt.x);
        lt.M_x = max(lt.x);
        lt.S_x = length(lt.x);
        lt.m_y = min(lt.y);
        lt.M_y = max(lt.y);
        lt.S_y = length(lt.y);
    end
    b1 = access_lt(log(la), log(nu), lt);
