function [I,W] = my_scol2im(Z,patchSize,sz,jmp)

patchSizeHat = patchSize*jmp - (jmp-1);

I = zeros(sz);
W = zeros(sz);

for idx=1:size(Z,2)
    patch = reshape(Z(:,idx),patchSize,patchSize);
    
    [i,j] = ind2sub(size(I)-patchSizeHat+1,idx);
    I(i:jmp:i+patchSizeHat-1,j:jmp:j+patchSizeHat-1) = I(i:jmp:i+patchSizeHat-1,j:jmp:j+patchSizeHat-1) + patch;
    W(i:jmp:i+patchSizeHat-1,j:jmp:j+patchSizeHat-1) = W(i:jmp:i+patchSizeHat-1,j:jmp:j+patchSizeHat-1) + ones(patchSize);
end

I = I./W;


