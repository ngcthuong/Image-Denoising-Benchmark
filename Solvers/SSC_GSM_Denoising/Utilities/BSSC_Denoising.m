function [im_out PSNR SSIM ]   =  BSSC_Denoising( par )
time0         =   clock;
nim           =   par.nim;
ori_im        =   par.I;
b             =   par.win;
[h  w ch]     =   size(nim);

N        =   h-b+1;
M        =   w-b+1;
r        =   [1:N];
c        =   [1:M]; 
disp(sprintf('PSNR of the noisy image = %f \n', csnr(nim, ori_im, 0, 0) ));

im_out      =   nim;
lamada      =   par.w;
nsig        =   par.nSig;

for iter = 1 : par.K        
    
    im_out  =   im_out + lamada*(nim - im_out);
    dif     =   im_out-nim;
    vd      =   nsig^2-(mean(mean(dif.^2)));
        
    if (iter==1)
        par.nSig  = sqrt(abs(vd));            
    else
        par.nSig  = sqrt(abs(vd))*par.lamada;
    end    
    
    if (mod(iter,6)==0) || (iter==1)
        [blk_arr,wei_arr]    =   Block_matching( im_out, par);
    end
    X         =   Im2Patch( im_out, par );
    
    Ys        =   zeros( size(X) );        
    W         =   zeros( size(X) );
    L         =   size(blk_arr,2);
    for  i  =  1 : L
        B          =   X(:, blk_arr(:, i));
        wei        =   repmat(wei_arr(:, i)',size(B, 1), 1);
        mB         =   repmat(sum(wei.*B, 2), 1, size(B, 2));
        %mB         =   repmat(mean( B, 2 ), 1, size(B, 2));
        B          =   B-mB;        
        %[Ys(:, blk_arr(:,i)), W(:, blk_arr(:,i))]   =   Low_rank_SSC( double(B), par.c1, par.nSig, mB );
        [Ys(:, blk_arr(:,i)), W(:, blk_arr(:,i))]   =   Bayes_SSC( double(B), par.c1, par.c2, par.nSig, mB );
    end

    im_out   =  zeros(h,w);
    im_wei   =  zeros(h,w);
    k        =   0;
    for i  = 1:b
        for j  = 1:b
            k    =  k+1;
            im_out(r-1+i,c-1+j)  =  im_out(r-1+i,c-1+j) + reshape( Ys(k,:)', [N M]);
            im_wei(r-1+i,c-1+j)  =  im_wei(r-1+i,c-1+j) + reshape( W(k,:)', [N M]);
        end
    end
    im_out  =  im_out./(im_wei+eps);
    
    if isfield(par,'I')
        PSNR      =  csnr( im_out, par.I, 0, 0 );
        SSIM      =  cal_ssim( im_out, par.I, 0, 0 );
        imwrite( im_out./255, 'Results\tmp.tif' );
    end
    
    fprintf( 'Iteration %d : nSig = %2.2f, PSNR = %2.2f, SSIM = %2.4f\n', iter, par.nSig, PSNR, SSIM );
end
if isfield(par,'I')
   PSNR      =  csnr( im_out, par.I, 0, 0 );
   SSIM      =  cal_ssim( im_out, par.I, 0, 0 );
end
disp(sprintf('Total elapsed time = %f min\n', (etime(clock,time0)/60) ));
return;



%------------------------------------------------------------------
% Re-weighted SV Thresholding
% Sigma = argmin || Y-U*Sigma*V' ||^2 + tau * || Sigma ||_*
%------------------------------------------------------------------
function  [X W]   =   Low_rank_SSC( Y, c1, nsig, m )
[U0,Sigma0,V0]    =   svd(full(Y),'econ');
Sigma0            =   diag(Sigma0);

S                 =   max( Sigma0.^2/size(Y, 2) - nsig^2, 0 );
thr               =   c1*nsig^2./ ( sqrt(S) + eps );
S                 =   soft(Sigma0, thr);
r                 =   sum( S>0 );

U                 =   U0(:,1:r);
V                 =   V0(:,1:r);
X                 =   U*diag(S(1:r))*V';

if r==size(Y,1)
    wei           =   1/size(Y,1);
else
    wei           =   (size(Y,1)-r)/size(Y,1);
end
W                 =   wei*ones( size(X) );
X                 =   (X + m)*wei;
return;


function  [X W]   =   Bayes_SSC( Y, c1, c2, nsig, mB )
[n, m]             =   size(Y);
[U,Lambda0,V0]     =   svd(full(Y),'econ');
theta0             =   sqrt( diag(Lambda0).^2/m );  % The initial \theta computed from noisy observations
B0                 =   sqrt(m)*V0';   % Compute B^(0)
A0                 =   Lambda0*V0';   % Compute A^(0)

% Compute \sigma^(0)
theta             =    sqrt( max( theta0.^2 - nsig^2, 0 ) );  % Compute theta^(0) using maximual likelihood method, used for computing \tau only.
tau               =    c1*nsig^2./ ( theta + eps );         
w                 =    diag( B0*B0' );
a                 =    diag( B0*A0' );
theta             =    soft(a, tau)./(w+eps);
t1                =    1./(theta.^2 + c2*nsig^2+eps);
B                 =    diag(t1)*diag(theta)'*A0;

% for  i = 1 : 1
%     w          =    diag( B*B' );
%     a          =    diag( B*A0' );
%     tau        =    c1*nsig^2./(theta+eps);
%     theta      =    soft(a, tau)./(w+eps);
%     t1         =    1./(theta.^2 + c2*nsig^2+eps);
%     B          =    diag(t1)*diag(theta)'*A0;    
% end

X                 =   U*diag( theta )*B;
 
% r                 =   sum( abs( diag(theta)*B )>0, 1);
% wei               =   min((size(Y,1)-r+1)/size(Y,1), 1);
% W                 =   repmat(wei, size(X,1), 1);

r                 =   sum( theta>0 );
if r==size(Y,1)
    wei           =   1/size(Y,1);
else
    wei           =   (size(Y,1)-r)/size(Y,1);
end
W                 =   wei*ones( size(X) );

X                 =   (X + mB).*W;
return;





% function  [X W]   =   Bayes_SSC( Y, c1, c2, nsig, m )
% [U0,Sigma0,V0]    =   svd(full(Y),'econ');
% Sigma             =   diag(Sigma0);
% 
% S                 =   max( Sigma.^2/size(Y, 2) - nsig^2, 0 );  % \sigma_i
% thr               =   c1*nsig^2./ ( sqrt(S) + eps );          
% S                 =   soft(Sigma, thr);
% 
% A0                =   Sigma0*V0';
% n                 =   size(V0, 1);
% Thr               =   repmat( c2*nsig^2./(S+eps),  1, n);
% A                 =   soft( A0, Thr );
% 
% 
% V                 =   (diag(1./(S+eps))*A);
% VtV               =   V*V';
% d                 =   diag(VtV);
% VtA               =   V*A0';
% a                 =   diag(VtA);
% 
% %S                 =   max( S.^2/size(Y, 2) - nsig^2, 0 );  % \sigma_i
% thr               =   c1*nsig^2./ (S + eps );   
% 
% S                 =   soft( a, thr )./(d+eps);
% 
% Thr               =   repmat( c2*nsig^2./(S+eps),  1, n);
% A                 =   soft( A0, Thr );
% 
% 
% 
% X                 =   U0*A;
% r                 =   sum( abs(A)>0, 1);
% 
% % r                 =   sum( S>0 );
% % U                 =   U0(:,1:r);
% % V                 =   V0(:,1:r);
% % X                 =   U*diag(S(1:r))*V';
% 
% wei               =   min((size(Y,1)-r+1)/size(Y,1), 1);
% W                 =   repmat(wei, size(X,1), 1);
% X                 =   (X + m).*W;
% return;
% 
