function y = solveInv_SR(x,imsize,lambda,weights,betas,patchSize,filters,BTB,tflag)

x = reshape(x,imsize);

% compute BTBX
y = lambda*BTB*x(:);

for i=1:length(filters)
    filtered = imfilter(x,filters{i});
    comp = imfilter(filtered,filters{i});
    y = y + weights(i)*betas(i)*patchSize^2*comp(:);
end


