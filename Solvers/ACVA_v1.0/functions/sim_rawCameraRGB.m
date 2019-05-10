function [R,G1,G2,B]=sim_rawCameraRGB(y0)
% From R,G,B to R, G1, G2, B

y_R=y0(:,:,1);
y_G=y0(:,:,2);
y_B=y0(:,:,3);
[m,n]=size(y_R);
m_even=2:2:m;
m_odd=1:2:m;
n_even=2:2:n;
n_odd=1:2:n;

mi=0;
for i=m_odd
    mi=mi+1;
    ni=0;
    for j=n_odd
        ni=ni+1;
        R(mi,ni)=y_R(i,j);
    end
   
end

mi=0;

for i=m_odd
    mi=mi+1;
    ni=0;
    for j=n_even
        ni=ni+1;
        G1(mi,ni)=y_G(i,j);
    end
   
end

mi=0;

for i=m_even
    mi=mi+1;
    ni=0;
    for j=n_odd
        ni=ni+1;
        G2(mi,ni)=y_G(i,j);
    end
   
end

mi=0;

for i=m_even
    mi=mi+1;
    ni=0;
    for j=n_even
        ni=ni+1;
        B(mi,ni)=y_B(i,j);
    end
   
end