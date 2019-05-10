function Z = my_im2col(I,patchSize,jmp)

patchSizeHat = patchSize*jmp - (jmp-1);

Z = zeros(patchSize^2,prod(size(I)-patchSizeHat+1));

for j=1:size(I,2)-patchSizeHat+1
    for i=1:size(I,1)-patchSizeHat+1
        patch = I(i:jmp:i+patchSizeHat-1,j:jmp:j+patchSizeHat-1);
        Z(:,sub2ind(size(I)-patchSizeHat+1,i,j)) = patch(:);
    end
end


