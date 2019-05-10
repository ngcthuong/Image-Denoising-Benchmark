function [cleanI,psnr] = sr(noiseI,cleanI,I,models,range,lambda,jmp,filters,weights,B,pad_sz)
% noiseI  - degraded image
% cleanI  - initialization of clean image
% I       - clean image (for computing PSNR)
% models  - cell of GMM models for the different scales
% range   - beta values
% lambda  - weight of degraded image (lambda = patchSize^2/noiseSD^2)
% jmp     - downsampling factors
% filters - blur filters applied before downsampling
% weights - the weights for the different scales
% B       - blur kernel used for degradation
% pad_sz  - size of padding (to avoid boundary issues)

% params
patchSize = 8;

BTY = B'*noiseI(:);
BTB = B'*B;

psnr = [];
for beta = range
    betas = zeros(1,length(filters));
    for i=1:length(betas)
        betas(i) = beta / norm(filters{i},'fro')^2;
    end
    
    % cell of all filters
    comp = cell(1,length(filters));

    parfor i=1:length(filters)
		% denoise patches

        % filter
		filtered = imfilter(cleanI,filters{i},'symmetric');

		% extract patches
		Z = my_im2col(filtered,patchSize,jmp(i));

		% remove mean
		meanZ = mean(Z);
		Z = bsxfun(@minus,Z,meanZ);
        
        SigmaNoise = 1/betas(i)*eye(patchSize^2);

		% calculate assignment probabilities for each mixture component for all patches
		PYZ = zeros(models{i}.nmodels,size(Z,2));
		for k=1:models{i}.nmodels
			PYZ(k,:) = log(models{i}.mixweights(k)) + loggausspdf2(Z,models{i}.covs(:,:,k) + SigmaNoise);
		end

		% find the most likely component for each patch
		[~,ks] = max(PYZ);

		% and now perform weiner filtering
		cleanZ = zeros(size(Z));
		for k=1:models{i}.nmodels
			inds = find(ks==k);
			cleanZ(:,inds) = (models{i}.covs(:,:,k) + SigmaNoise) \ (models{i}.covs(:,:,k) * Z(:,inds));
		end

		% add back mean
		cleanZ = bsxfun(@plus,cleanZ,meanZ);

		% upsample and estimate clean image
		upsampled = my_scol2im(cleanZ,patchSize,size(cleanI),jmp(i));

		% backward filter
		comp{i} = imfilter(upsampled,filters{i},'symmetric','conv');
    end

    % conjugate gradient
    b = lambda*BTY;
    for i=1:length(filters)
        b = b + weights(i)*betas(i)*patchSize^2*comp{i}(:);
    end
    
    aux_fun = @(x,tflag) solveInv_SR(x,size(cleanI),lambda,weights,betas,patchSize,filters,BTB,tflag);
    [tmp,~] = bicg(aux_fun,b,1e-6,10000);
    cleanI = reshape(tmp,size(cleanI));

    % remove pad
    I_ = I(1+pad_sz:end-pad_sz,1+pad_sz:end-pad_sz);
    cleanI_ = cleanI((1+pad_sz):(pad_sz+size(I_,1)),(1+pad_sz):(pad_sz+size(I_,2)));

    % psnr
    psnr = [psnr 20*log10(1/std2(cleanI_-I_))];

%     figure(1); plot(psnr); title('PSNR');
%     drawnow;
end




