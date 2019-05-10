function est13 = progressiveTHR(X,level,scal)


% The progressive denoising 


if  isempty(X)
    est13=0;
else
% MPSVD denoising
    [M,L]=size(X);
    gamma=M/L;
    thresh=scal*level*level*(1+sqrt(gamma))*(1+sqrt(gamma));%

    X0=X';
    mx=mean(X0(:));

    X0=X0-mx;
    if L>30
        [V,D] = eig(X0'*X0/(L-1));
    else
        [V,D] = eig(X0'*X0/(L));
    end

    Dc=D-eye(size(D))*thresh;

    D_m=max(Dc,zeros(size(D)));
    D_diag=diag(D_m);
    D_k=sum(D_diag>0);

    nb_axis=max(D_k,1);

    V = V(:,(end-nb_axis+1):end);%
    S = sqrt(D((end-nb_axis+1):end,(end-nb_axis+1):end));%
    U = X0*V*S^(-1);  

% LMMSE denoising
    P = U*S;
    lev = level*level;
    px = P.^2;


    if size(px,1)>3

        px = [px(1,:);(px(2:end-1,:)+px(1:end-2,:)+px(3:end,:))/3;px(end,:)];

    else

        px = repmat(mean(px,2),1,size(px,2));   
    end


    px0   =   max(0, px-lev);
    weight  =   px0./px;         
    weight(isnan(weight))=0;
    est13 = (P.*weight)*V'+mx;


    est13 = est13';
    est13 = est13.*(est13>0); 
    
end 
