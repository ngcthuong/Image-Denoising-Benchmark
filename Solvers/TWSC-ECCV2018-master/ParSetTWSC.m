function Par = ParSetTWSC(nSig)

Par.image = 1;
Par.nSig = nSig/255;
Par.innerIter = 2;
Par.win = 30;
Par.lambda1 = 0;
Par.ps = 8;
Par.outerIter = 10;
Par.step = 3;
Par.nlspini = 90;
Par.nlspgap = 10;
if 0 < nSig <= 20
    Par.outerIter = 8;
    Par.delta = .07;
    Par.nlspini = 70;
    Par.lambda2 = .9;
elseif 20 < nSig <= 30
    Par.delta = .06;
    Par.lambda2 = .76;
elseif 30 < nSig <= 40
    Par.delta   = .07;
    Par.lambda2 = .78;
elseif 40 < nSig <= 60
    Par.nlspini = 120;
    Par.nlspgap = 15;
    Par.delta   = .05;
    Par.lambda2 = .72;
elseif 60 < nSig <= 80
    Par.ps = 9;
    Par.outerIter = 14;
    Par.step = 4;
    Par.nlspini = 140;
    Par.delta = .05;
    Par.lambda2 = .68; % .66
else
    disp('Please tune the above parameters by yourself, thanks!');
end
