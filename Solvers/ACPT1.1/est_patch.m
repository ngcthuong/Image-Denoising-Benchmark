function sigma = est_patch(rima,Patch_width)

data1=rima;
d=Patch_width;
stride=1;
X=image2cols(data1,d,stride);  

[M,L]=size(X);
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
sigma=sqrt(sigmaest1(ind));

