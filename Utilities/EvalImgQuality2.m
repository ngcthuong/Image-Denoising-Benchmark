% Author: Shunsuke Ono (ono@sp.ce.titech.ac.jp)

function [ePSNR, eSSIM] = EvalImgQuality2(ImgRec, ImgOrg)
ePSNR = calPSNR(ImgRec, ImgOrg);
eSSIM = cal_ssim( ImgRec, ImgOrg, 0, 0 );


