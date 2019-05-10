function sigma = est_x(X)
% Noise level estimation
% The algorithm is described in the following article:
% "Detail-preserving Image Denoising via Adaptive Clustering and Progressive PCA Thresholding"
% Wenzhao Zhao, Yisong Lv, Qiegen Liu*, and Binjie Qin*. IEEE Access.

if  isempty(X)
    sigma=0;
else
    [M,L]=size(X);
    if M>L
        X=X';
        [M,L]=size(X);
    end
    mx=mean(X, 1);
    X=X-repmat(mx,[M, 1]);

    [~,Eg]=eig(X*X'/(L));
    Eg=diag(Eg);

    sigmaest1=zeros(M-1,1);
    sigmaest2=zeros(M-1,1);
    for i=M:-1:2

        sigmaest1(i-1)=sum(Eg(1:i,:))/i;
        sigmaest2(i-1)=(Eg(i)-Eg(2))/4/(sqrt((M)/L)); 

    end
    [~,ind]=min(abs(sigmaest1-sigmaest2));
    %
    sigma=sqrt(abs(sigmaest1(ind)));%abs
end
