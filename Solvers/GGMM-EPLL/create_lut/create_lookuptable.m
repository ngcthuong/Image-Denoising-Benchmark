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

clear all
close all

cd(fileparts(mfilename('fullpath')));

desc = 'x=log(lambda), y=log(nu)';

load('gamma0_lookup_table.mat')
load('a1_lookup_table.mat')
load('a2_lookup_table.mat')
load('b1_lookup_table.mat')
load('b2_lookup_table.mat')
load('h_lookup_table.mat')

lut.log_ga0 = log_ga0;
lut.a1  = a1;
lut.a2  = a2;
lut.b1  = b1;
lut.b2  = b2;
lut.h   = h;
lut.x   = log_la;
lut.y   = log_nu;
lut.M_x = max(lut.x);
lut.m_x = min(lut.x);
lut.M_y = max(lut.y);
lut.m_y = min(lut.y);
lut.S_x = length(lut.x);
lut.S_y = length(lut.y);

save('params_lookup_table.mat', 'desc', 'lut');