function par  =  ParSetGSRC( nSig,I )

par.I                =      double(I);

par.nSig            =       nSig;

par.Iter             =       20;

par.Thr              =       0.0001;

par.step             =        4;

if nSig <=30

par.c                =        0.2*2.0*sqrt(2);  
 
par.w                =        0.18;

par.lamada           =        0.67;

par.win              =         7;

par.nblk             =         60;

par.Thr              =       0.0001;

elseif nSig <=40

par.c                =        0.3*2.0*sqrt(2);  

par.w                =        0.22;

par.lamada           =        0.67;

par.win              =         7;

par.nblk             =         60;

par.Thr              =        0.00008;

elseif nSig <=50

par.c                =        0.3*2.0*sqrt(2);  

par.w                =        0.22;

par.lamada           =        0.67;

par.win              =         7;

par.nblk             =         60;

par.Thr              =       0.0001;

elseif nSig <=75

par.c                =        0.3*2.0*sqrt(2);  

par.w                =        0.22;

par.lamada           =        0.67;

par.win              =         8;

par.nblk             =         80;

par.Thr              =       0.0001;

else
    
par.c                =        0.3*2.0*sqrt(2);  

par.w                =        0.22;

par.lamada           =        0.67;

par.win              =         9;

par.nblk             =         90;

par.Thr              =       0.0001;

end

end

