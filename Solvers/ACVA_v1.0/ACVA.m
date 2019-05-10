%% Source code for submission to IEEE TIP 2018
%% Author : Wenzhao Zhao
%% Matlab version : R2013a or R2015b

function [RES,RES0,NUM]=ACVA(rima,wid,step,level,seed)
if (exist('level') ~= 1)
    level='n';
end
if (exist('seed') ~= 1)
    seed='n';
end

[M,N] = size(rima);
RES0 = zeros(M,N);
NUM = zeros(M,N);

 range_x = 1:step:(M-wid+1);
 range_y = 1:step:(N-wid+1);
  if (range_x(end)~=(M-wid+1))
    range_x = [range_x (M-wid+1)];
  end
  if (range_y(end)~=(N-wid+1))
    range_y = [range_y (N-wid+1)];
  end
  
totalnum = length(range_y)*length(range_x);
countnum = 0;
ProgressBar(countnum,totalnum);
for  x = range_x%
    
    for y = range_y%

        [~,cres0,num0] = CMPLRA_bands(rima(x:x+wid-1,y:y+wid-1),level,seed);% denoising in the sliding window 
        countnum = countnum+1;
        ProgressBar(countnum,totalnum);
        %disp(countnum/totalnum);
        
        RES0(x:x+wid-1,y:y+wid-1) = RES0(x:x+wid-1,y:y+wid-1) + cres0;
        NUM(x:x+wid-1,y:y+wid-1) = NUM(x:x+wid-1,y:y+wid-1) + num0;
        
    end
    
end

RES=RES0./NUM;

function [Y1,Y10,num0]= CMPLRA_bands(rima,level,seed)

%% source code for submission to IEEE TIP 
%% author : Wenzhao Zhao


data = rima;

data1 = data;%
[m,n] = size(data1);

Patch_width = 8;
d = Patch_width;
s = size(rima);
stride = 1;
X = image2cols(data1, d, stride);   %extracting image patches
[M,N] = size(X);

if level == 'n'
   level = est_x(X);
end

Xc = X;%
%% the first stage clustering

k44 = max(floor(m*n/(256*256)),4);

if (seed~= 'n')
label4 = litekmeans_k(Xc, k44,seed); 
else
label4 = litekmeans_k(Xc, k44); 
end


indexes4 = cell(k44,1);
for i = 1:k44
    indexes4{i} = find(label4 == i); 
end 
%% the second stage clustering
labels = cell(k44,1);
knm = zeros(k44,1);
for i = 1:k44
    knm(i) = max(floor(length(indexes4{i})/d/d),1);

    if (seed ~= 'n')
        labels{i} = litekmeans_k(Xc(:,indexes4{i}), knm(i),seed); 
    else
        labels{i} = litekmeans_k(Xc(:,indexes4{i}), knm(i)); 
    end
end   

k = sum(knm);

indexes = cell(k,1);
ijindex = 0;
for i = 1:k44
    for j = 1:knm(i)
        ijindex = ijindex+1;

        tempindex = (labels{i} == j);
        indexes{ijindex} = (indexes4{i}(tempindex));
    end 
end

mean_cluster = zeros(size(Xc,1),k); 

for i = 1:k
    size_cluster(i) = size((indexes{i}),1);
    mean_cluster(:,i) = mean(Xc(:,indexes{i}),2);
end

 
%% merging small clusters 
thresh = level*level*0.25;%
 
index_new = combination_size(mean_cluster,size_cluster,indexes,k,thresh,X);
k1 = size(index_new,2);
mean_cluster_n = zeros(size(Xc,1),k1); 
for i = 1:k1
    size_cluster_n(i) = size((index_new{i}),1);
    mean_cluster_n(:,i) = mean(Xc(:,index_new{i}),2);
end

k2 = 0;

while k1 ~= k2
    k1 = size(index_new,2);
 
    index_new = combination_size(mean_cluster_n,size_cluster_n,index_new,k1,thresh,X);
    k2 = size(index_new,2);
    mean_cluster_n = zeros(size(Xc,1),k2); 
    for i = 1:k2
        size_cluster_n(i) = size((index_new{i}),1);
        mean_cluster_n(:,i) = mean(Xc(:,index_new{i}),2);
    end
%     
end

%% the two-step denoising
denoi_stacked1 = zeros(size(Xc));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
for i = 1:k2

    input_x=X(:,index_new{i});

    scal1 = 1.1;%
    est1 = MPLRA_LMMSE_LPA(input_x,level,scal1);

% % % % % % % % % % % 
    if ~isempty(index_new{i})
        denoi_stacked1(:,index_new{i}) = est1;
    end
    
end


[Y1,Y10,num0] = columns2im(denoi_stacked1, s, stride);



function label = litekmeans_k(X, k,seed)
% Perform k-means clustering.
%   X: d x n data matrix
%   k: number of seeds
% Written by Michael Chen (sth4nth@gmail.com).
% Copyright (c) 2009, Michael Chen
% All rights reserved.

% Modified by Wenzhao Zhao

i = 0;
n = size(X,2);

if k == 1
    label = ones(1,n);
	if verLessThan('matlab', '8.1') == 0 % transposes depending on MATLAB version
		label = label';
	end
else

	last = 0;
if (exist('seed') == 1)
	rand('state',seed);
end
	label = ceil(k*rand(1,n));  % random initialization
	while any(label ~= last)
		[~,~,label] = unique(label);   % remove empty clusters
		E = sparse(1:n,label,1,n,k,n);  % transform label into indicator matrix
		center = X*(E*spdiags(1./sum(E,1)',0,k,k));    % compute center of each cluster
		last = label;
		[~,label] = max(bsxfun(@minus,center'*X,0.5*sum(center.^2,1)')); % assign samples to the nearest centers
		if verLessThan('matlab', '8.1') == 0 % transposes depending on MATLAB version
			label = label';
		end
		i = i+1;
 
	   if i > 200
		   break;
	   end
	end

end




function index_new = combination_size(mean_cluster,size_cluster,indexes,k,thresh,X)
%%  merging clusters
M = size(X,1);
if k == 1
    index_new = indexes;
else
    min2 = zeros(1,k);
    min_i = zeros(1,k);
    for i = 1:k


        minratio = 1-min(size_cluster)/size(X,2);

        if i ~= k
            dif_m = sum((repmat(mean_cluster(:,i),1,k-i)-mean_cluster(:,i+1:k)).^2)/M;

            minsize = min(size_cluster(i),size_cluster(i+1:k));
 
            siz_m = (minsize>0)-1.*(0.3*(minsize>200));%               
 
            dif_m = dif_m./siz_m;
            [dif,min_i11] = min(dif_m);
            min_i(i) = min_i11+i;
        else
            min_i(i) = k-1; 
            dif = sum((mean_cluster(:,k)-mean_cluster(:,k-1)).^2)/M;
 
            minsize = min(size_cluster(k),size_cluster(k-1));
 

            ssm = (minsize>0)-1.*(0.3*(minsize>200));% 
 
            dif = dif/ssm;
            
 
        end
        min2(1,i) = dif;
    end

    ind_re = 1:k;
    for i = 1:k
        if min2(i) < thresh 
            ind_re(min_i(i)) = min(i,ind_re(i));
        end


    end
    te = 0;
    index_new = {};
    for i = 1:k
        II = find(ind_re == i);
        if size(II)
 
            te = te+1;
            in=[];
            for tt = 1:length(II)
               in = [in ; indexes{II(tt)}];
            end
            index_new{te} = in;
        end
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function est13 = MPLRA_LMMSE_LPA(X,level,scal)
%% Matrix denoising
X(isinf(X)) = 0;
X(isnan(X)) = 0;
% The two-step denoising 

if size(X,1) == 0
    warning('off')
end
if  isempty(X)
    est13=0;
else

    [M,L] = size(X);
    gamma = M/L;
    thresh = scal*level*level*(1+sqrt(gamma))*(1+sqrt(gamma));% Threshold

    X0 = X';
    mx = mean(X0(:));

    X0 = X0-mx;
 
    [V,D] = eig(X0'*X0/(L));
 
    Dc = D-eye(size(D))*thresh;
    D_m = max(Dc,zeros(size(D)));   
    
    D_diag = diag(D_m);
    D_k = sum(D_diag>0);
    
    nb_axis = max(D_k,1);
    
    S = sqrt(D((end-nb_axis+1):end,(end-nb_axis+1):end));
 
    V = V(:,(end-nb_axis+1):end);

    U = X0*V*S^(-1);  

    P = U*S;
    lev = level*level;% 

	if level < 15
		Hw = [1 2 3];%  
	elseif level < 25;
		Hw = [1 2 3 4 5];% 
	else
		Hw =[3 5 7 10 15];%
	end
	gamma = 1;
	[~,px] = LPA_ICI_all(P',Hw,level,gamma); % LPA-ICI
	 
	px = px';
	  
	weight = subWiener(px,lev,0.7);% The suboptimal Wiener filter
    est13 = (P.*weight)*V'+mx;   
    est13 = est13'; 
    est13(isnan(est13)) = 0;
end 
if size(X,1) == 0
    warning('on')
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function weight = subWiener(px,sigma2,beta)
alp = 0:0.01:1;
alp2 = alp.^2;
mpx = px;
 
w = zeros(size(mpx));
for i = 1:size(w,1)
    
    for j = 1:size(w,2)
        b = sigma2/mpx(i,j);
        J = (1+b^2-2*b)./(1+alp2*b*b-2*alp*b)-beta*alp2;
        [~,Ind] = max(J);
        alpha = (Ind-1)*0.01;
        w(i,j) = max(1-alpha*b,0);
    end
    
end
weight = w;% 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     

function [P_est,px] = LPA_ICI_all(P,Hw,level,gam)
% Hw: window sizes

N = size(P,2);
P_est = zeros(size(P));
 
px = zeros(size(P));

for nm = 1:size(P,1)
    p1d = P(nm,:);
    for si = 1:N

        [p_est,~, pxn] = LPA_ICI(p1d,si,N,Hw,level,gam);
        P_est(nm,si) = p_est;
        px(nm,si) = pxn;
        
    end

end

function [p_est,h_opt,px] = LPA_ICI(P,si,N,Hw,level,gamma)
% P is of size 1*N

pxi = zeros(1,length(Hw));
Ch = zeros(1,length(Hw));
std = zeros(length(Hw),1);

hi = 0;
for h=Hw 
      
    hi = hi+1;
    dx = max(1-si,-h):min(h,N-si); 
    ddx = si+max(1-si,-h):si+min(h,N-si); 
    z = P(:,ddx);
    X = dx;

    b0 = length(X);
    bz = sum(z);
    C = bz/b0;
    Ch(:,hi) = C;

    pxi(hi) = mean(z.*z);

    std(hi) = level/sqrt(b0);%
       
end

Band = Ch(1,:)';
[p_est, IND, ~] = apply_ICI(Band,std,gamma);
h_opt = Hw(IND);
px = pxi(IND);
 
function [est, ind, std_ad] = apply_ICI(band,std,gamma)

est = band(1);
std_ad = std(1);
intv(1) = 1;
ind = 1;

D = gamma*std;
U = band+D;
L = band-D;
len = length(band);

for i = 2:len
    
    L(i) = max(L(i-1),L(i));         
    U(i) = min(U(i-1),U(i));
    intv(i) = (U(i) >= L(i));                  
    intv(i) = intv(i) & intv(i-1);
    
    est = intv(i)*band(i)+~intv(i)*est;  
    ind = intv(i)*i+~intv(i)*ind;  
    std_ad = intv(i)*std(i)+~intv(i)*std_ad; 
    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function res = image2cols(im, pSz, stride)

  range_y = 1:stride:(size(im,1)-pSz+1);
  range_x = 1:stride:(size(im,2)-pSz+1);
  if (range_y(end) ~= (size(im,1)-pSz+1))
    range_y = [range_y(size(im,1)-pSz+1)];
  end
  if (range_x(end) ~= (size(im,2)-pSz+1))
    range_x = [range_x(size(im,2)-pSz+1)];
  end
  sz = length(range_y)*length(range_x);

  res = zeros(pSz^2, sz);

  idx = 0;
  for y = range_y
    for x = range_x
      p = im(y:y+pSz-1,x:x+pSz-1);
      idx = idx + 1;
      res(:,idx) = p(:);
    end
  end

function [res,res0,w] = columns2im(cols, im_sz, stride)

  pSz = sqrt(size(cols,1));
  res0 = zeros(im_sz(1), im_sz(2));
  w = zeros(im_sz(1), im_sz(2));

  range_y = 1:stride:(im_sz(1)-pSz+1);
  range_x = 1:stride:(im_sz(2)-pSz+1);
  if (range_y(end)~=(im_sz(1)-pSz+1))
    range_y = [range_y(im_sz(1)-pSz+1)];
  end
  if (range_x(end)~=(im_sz(2)-pSz+1))
    range_x = [range_x(im_sz(2)-pSz+1)];
  end

  idx = 0;
  for y=range_y
    for x=range_x
      idx = idx + 1;
      p = reshape(cols(:,idx), [pSz pSz]);
      res0(y:y+pSz-1, x:x+pSz-1) = res0(y:y+pSz-1, x:x+pSz-1) + p;
      w(y:y+pSz-1, x:x+pSz-1) = w(y:y+pSz-1, x:x+pSz-1) + 1;
      
    end
  end
  res = res0./w;
