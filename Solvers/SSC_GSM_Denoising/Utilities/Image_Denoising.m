%----------------------------------------------------
% Blur: 9x9 Uniform kernel; AGWN std Varance = 2^0.5
% Data: May 20th, 2010
% Author: Weisheng Dong, wsdong@mail.xidian.edu.cn
%----------------------------------------------------
function  Image_Denoising( nSig, Out_dir, In_dir, idx )
pre           =   'BSSC_';
par.nSig      =   nSig;
levels        =   [5, 10,  15, 20, 50, 100];
K             =   [4,  6,  10, 11, 16, 16];
par.K         =   K(idx);

if nSig<=5
    par.win       =   6;
    par.nblk      =   43;
    par.c1        =   2.1*sqrt(2);    %3.2
    par.c2        =   par.c1/4;       %0.78
    
    par.lamada    =   0.62;  % .63
    par.w         =   0.32;  % .28
    par.hp        =   55;
    par.p         =   1;
elseif nSig<=15
    par.win       =   6;
    par.nblk      =   44;
    par.c1        =   3.2*sqrt(2);    %3.2    2.3
    par.c2        =   par.c1/4;       %0.78
    
    par.lamada    =   0.46;  % .63
    par.w         =   0.1;  % .32
    par.hp        =   70;
    par.p         =   1;
elseif nSig <= 20
    par.win       =   7;
    par.nblk      =   58;
    par.c1        =   3.1*sqrt(2);   % 2.9
    par.c2        =   par.c1/4;      % 0.7*sqrt(2);
    
    par.lamada    =   0.65;  %0.69
    par.w         =   0.23;
    par.hp        =   85;        
    par.p         =   0.67;
elseif nSig <= 30
    par.win       =   7;
    par.nblk      =   58;
    par.c1        =   3.1*sqrt(2);   % 2.9
    par.c2        =   par.c1/4;      % 0.7*sqrt(2);
    
    par.lamada    =   0.67;  %0.69
    par.w         =   0.23;
    par.hp        =   80;    
    par.p         =   0.67;
elseif nSig<=50
    par.win       =   8;   % 8
    par.nblk      =   76;  % 76
    par.c1        =   2.9*sqrt(2);   % 3.0
    par.c2        =   par.c1/4;
    
    par.lamada    =   0.67;
    par.w         =   0.23;    
    par.hp        =   90;
    par.p         =   0.67;
else
    par.win       =   9;   % 10
    par.nblk      =   90;  % 109
    par.c1        =   2.9*sqrt(2);   % 3.1
    par.c2        =   par.c1/4;   % 3.7
    
    par.lamada    =   0.67;
    par.w         =   0.23;    
    par.hp        =   95;
    par.p         =   0.67;
end
par.step      =   min(4, par.win-1);
par.sel       =   'random';
fpath         =   fullfile(In_dir, '*.tif');
im_dir        =   dir(fpath);
im_num        =   length(im_dir);
cnt           =   0;
sum_psnr      =   0;
sum_ssim      =   0;
time0         =   clock;
fn_txt        =   strcat( pre, 'PSNR_SSIM.txt' ); 
fd_txt        =   fopen( fullfile(Out_dir, fn_txt), 'wt');

for i = 1:im_num 
    
    par.I        =   double( imread(fullfile(In_dir, im_dir(i).name)) );
    par.nim      =   par.I + nSig*Gen_noise(In_dir, im_dir, i);
    imwrite(par.nim./255,'tmp.tif');
    
    [im PSNR SSIM]   =   BSSC_Denoising2( par );
    sum_psnr    =  sum_psnr + PSNR;
    sum_ssim    =  sum_ssim + SSIM;    
    
    fname            =   strcat(pre, im_dir(i).name);
    imwrite(im./255, fullfile(Out_dir, fname));
    disp( sprintf('%s: PSNR = %3.2f  SSIM = %f\n', im_dir(i).name, PSNR, SSIM) );
    fprintf(fd_txt, '%s :  PSNR = %2.2f  SSIM = %2.4f\n', im_dir(i).name, PSNR, SSIM);
    cnt   =  cnt + 1;
end
fprintf(fd_txt, '\n\nAverage :  PSNR = %2.2f  SSIM = %2.4f\n', sum_psnr/cnt, sum_ssim/cnt);
fclose(fd_txt);
disp(sprintf('Total elapsed time = %f min\n', (etime(clock,time0)/60) ));
return;


function nim  =  Gen_noise( In_dir, im_dir, i )
randn('seed',0);
for ii=1:i
    im        =   imread(fullfile(In_dir, im_dir(ii).name));
    nim       =   randn(size(im));
end
return;

