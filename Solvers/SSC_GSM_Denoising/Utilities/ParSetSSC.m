function  par  =  ParSetSSC( nSig, K )
par.nSig      =   nSig;
%K             =   [3,  3,  3,  3, 4, 5];
par.K         =   K;
par.sigma     =   1.7;

if nSig<=5
    par.win       =   6;
    par.nblk      =   43;
    par.c1        =   2.0*sqrt(2);
    par.c2        =   par.c1/4;
    
    par.lamada    =   0.62;
    par.w         =   0.32;
    par.hp        =   60;
elseif nSig<=15
    par.win       =   6;
    par.nblk      =   44;
    par.c1        =   1.9*sqrt(2);
    par.c2        =   par.c1/4; 
    
    par.lamada    =   0.63; 
    par.w         =   0.30; 
    par.hp        =   70;
elseif nSig <= 20
    par.win       =   7;
    par.nblk      =   50;   
    par.c1        =   1.5*sqrt(2);  
    par.c2        =   par.c1/4; 
    
    par.lamada    =   0.67;
    par.w         =   0.21; 
    par.hp        =   80;        
elseif nSig <= 30
    par.win       =   7;
    par.nblk      =   58;
    par.c1        =   2.0*sqrt(2);  
    par.c2        =   par.c1/4;  
    
    par.lamada    =   0.67;
    par.w         =   0.22;
    par.hp        =   80;    
elseif nSig<=50
    par.win       =   7;
    par.nblk      =   58;
    par.c1        =   1.8*sqrt(2);  
    par.c2        =   par.c1/4;  
    
    par.lamada    =   0.67;
    par.w         =   0.22;
    par.hp        =   80;    
else
    par.win       =   9; 
    par.nblk      =   90;
    par.c1        =   2.1*sqrt(2);
    par.c2        =   par.c1/4; 
    
    par.lamada    =   0.67;
    par.w         =   0.23;    
    par.hp        =   95;
end
par.step      =   min(4, par.win-1);
par.cls_num   =   50;
