function Y1= ACPT(rima,level)
% 
% % % source code for submission to IEEE Access : Detail-preserving Image Denoising via Adaptive Clustering and Progressive PCA Thresholding
% % 
% 
% % noise level estimation
if (exist('level') ~= 1),
level=noiselevel(rima);

end

% % image denoising with ACPT
data=rima;
Mi=min(data(:));
data1=data-Mi;
[m,n]=size(data1);

Patch_width=8;%
d=Patch_width;
s=size(rima);
stride=1;
X = image2cols(data1, d, stride);   %extracting image patches

% % the first stage clustering
[M,N]=size(X);
k44=max(floor(m*n/(256*256)),4);
label4 = litekmeans_m(X, k44);  


indexes4=cell(k44,1);
for i = 1:k44
    indexes4{i} = find(label4 == i); 
end 
% % the second stage clustering
labels=cell(k44,1);
knm=zeros(k44,1);
for i = 1:k44
    knm(i)=max(floor(length(indexes4{i})/d/d),1);
    labels{i} = litekmeans_m(X(1:end,indexes4{i}), knm(i)); 

end   

k=sum(knm);
disp(['Overclustering: cluster number : ' num2str(k) ]);
indexes=cell(k,1);
ijindex=0;
for i=1:k44
    for j = 1:knm(i)
        ijindex=ijindex+1;
        tempindex=find(labels{i}==j);
        indexes{ijindex}=(indexes4{i}(tempindex));
    end 
end

mean_cluster = zeros(size(X,1),k); 

for i = 1:k
    mean_cluster(:,i) = mean(X(1:end,indexes{i})')';
end

 disp('Iterative merging :');
% % merging small clusters 
thresh=level*level*0.25;

index_new=merging_thr(mean_cluster,indexes,k,thresh,X);
k1=size(index_new,2);
size_cluster_n = zeros(k1,1); 
mean_cluster_n = zeros(size(X,1),k1); 
for i = 1:k1
    
    size_cluster_n(i) = size((index_new{i}),1);
    mean_cluster_n(:,i) = mean(X(1:end,index_new{i})')';
end
k2=0;

while k1~=k2
    k1=size(index_new,2);
    disp(['cluster number : ' num2str(k1) ]);
    index_new=merging_thr(mean_cluster_n,index_new,k1,thresh,X);

    k2=size(index_new,2);
    size_cluster_n = zeros(k2,1); 
    mean_cluster_n = zeros(size(X,1),k2); 
    for i = 1:k2

        size_cluster_n(i) = size((index_new{i}),1);
        mean_cluster_n(:,i) = mean(X(1:end,index_new{i})')';
    end

end


disp('Progressive denoising ...');
% % the progressive denoising
denoi_stacked1 = zeros(size(X));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
for i = 1:k2


    input_x=X(1:end,index_new{i});
    
    scal1= 1.1;
    est1 = progressiveTHR(input_x,level,scal1);

    denoi_stacked1(:,index_new{i}) = est1;
    
end


Y1 = columns2im(denoi_stacked1, s, stride);
Y1=Y1+Mi;
