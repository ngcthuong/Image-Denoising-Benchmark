function  [par, model]  =  ParSetPGPD( nSig )

par.step = 3;       % the step of two neighbor patches
par.IteNum = 4;  % the iteration number
par.nSig      =   nSig/255;

if nSig <= 10
    load 'PG_GMM_6x6_win15_nlsp10_delta0.002_cls65.mat';
    par.c1 = 0.33*2*sqrt(2);
    par.delta = 0.1;
    par.eta=0.79;
elseif nSig<=20
    load 'PG_GMM_6x6_win15_nlsp10_delta0.002_cls65.mat';
    par.c1 = 0.29*2*sqrt(2);
    par.delta = 0.09;
    par.eta=0.73;
elseif nSig <=30
    load 'PG_GMM_7x7_win15_nlsp10_delta0.002_cls33.mat';
    par.c1 = 0.19*2*sqrt(2);
    par.delta = 0.08;
    par.eta=0.89;
elseif nSig<=40
    load 'PG_GMM_8x8_win15_nlsp10_delta0.002_cls33.mat';
    par.c1 = 0.15*2*sqrt(2);
    par.delta = 0.07;
    par.eta=0.98;
elseif nSig<=50
    load 'PG_GMM_8x8_win15_nlsp10_delta0.002_cls33.mat';
    par.c1 = 0.12*2*sqrt(2);
    par.delta = 0.06;
    par.eta=1.05;
elseif nSig<=75
    load '/PG_GMM_9x9_win15_nlsp10_delta0.002_cls33.mat';
    par.c1 = 0.09*2*sqrt(2);
    par.delta = 0.05;  
    par.eta=1.15;
else
    load '/PG_GMM_9x9_win15_nlsp10_delta0.002_cls33.mat';
    par.c1 = 0.06*2*sqrt(2);
    par.delta = 0.05;
    par.eta=1.30;
end
par.ps = ps;        % patch size
par.nlsp = nlsp;  % number of non-local patches
par.Win = win;   % size of window around the patch
% dictionary and regularization parameter
for i = 1:size(GMM_D,2)
    par.D(:,:,i) = reshape(single(GMM_D(:, i)), size(GMM_S,1), size(GMM_S,1));
end
par.S = single(GMM_S);
