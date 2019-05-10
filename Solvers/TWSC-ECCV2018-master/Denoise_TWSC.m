function [imRec, Par] = Denoise_TWSC(nim, Par)
% Function for Denoising TWSC 

Par.nim  = nim;

Par.PSNR = zeros(Par.outerIter, 1, 'double');
Par.SSIM = zeros(Par.outerIter, 1, 'double');
T512     = [];
T256     = [];

Par.nlsp = Par.nlspini;  % number of non-local patches
[imRec, Par]  =  TWSC_Sigma_AWGN(Par);

if size(Par.I,1) == 512
    T512 = [T512 etime(clock,time0)];
    fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );
elseif size(Par.I,1) ==256
    T256 = [T256 etime(clock,time0)];
    fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );
end
im_out(im_out>1)=1;
im_out(im_out<0)=0;


PSNR = Par.PSNR(end,:);
Par.mPSNR= mean(PSNR,2);
SSIM = Par.SSIM(end,:);
Par.mSSIM= mean(SSIM,2);
Par.mT512 = mean(T512);
Par.sT512 = std(T512);
Par.mT256 = mean(T256);
Par.sT256 = std(T256);

fprintf('The average PSNR = %2.4f, SSIM = %2.4f. \n', Par.mPSNR, Par.mSSIM);

%name = sprintf([write_MAT_dir method '_nSig' num2str(nSig) '_oIte' num2str(Par.outerIter) '_iIte' num2str(Par.innerIter) '_ps' num2str(Par.ps) '_step' num2str(Par.step) '_nlspini' num2str(Par.nlspini) '_nlspgap' num2str(Par.nlspgap) '_delta' num2str(Par.delta) '_l1' num2str(Par.lambda1) '_l2' num2str(Par.lambda2) '.mat']);
%save(name,'nSig','PSNR','SSIM','mPSNR','mSSIM','mT512','sT512','mT256','sT256');

end