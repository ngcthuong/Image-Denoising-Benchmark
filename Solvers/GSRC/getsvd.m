function [P, V]=getsvd(X)

%X: MxN matrix (M dimensions, N trials)
%Y: Y=P*X
%P: the transform matrix
%V: the variance vector

[M,N]=size(X);

mx   =  mean(X,2);
mx2  =  repmat(mx,1,N);
X    =  X-mx2;

CovX=X*X'/(N-1);

%CovX  =  CovX - diag(nsig^2*ones(size(CovX,1),1));

[P,V,~]=svd(CovX);

V=diag(V);

return;

