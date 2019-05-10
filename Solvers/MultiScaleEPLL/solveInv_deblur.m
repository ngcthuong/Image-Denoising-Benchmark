function y = solveInv_deblur(x,K,imsize,lambda,weights,betas,W,filters)

x = reshape(x,imsize);

% compute ATAX
tt = imfilter(x,K,'conv','same');
tt = tt(floor(size(K,1)/2)+1:end-floor(size(K,1)/2),floor(size(K,2)/2)+1:end-floor(size(K,2)/2));
y = lambda*imfilter(tt,rot90(rot90(K)),'conv','full');

for i=1:length(filters)
    filtered = imfilter(x,filters{i},'symmetric');
    comp = imfilter(filtered,filters{i},'symmetric','conv');
    y = y + weights(i)*betas(i)*W{i}.*comp;
end

y = y(:);


