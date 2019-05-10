function y = solveInv_denoise(x,imsize,lambda,weights,betas,patchSize,filters)

x = reshape(x,imsize);

y = lambda*x;
for i=1:length(filters)
    filtered = imfilter(x,filters{i},'symmetric');
    comp = imfilter(filtered,filters{i},'symmetric','conv');
    y = y + weights(i)*betas(i)*patchSize^2*comp;
end

y = y(:);


