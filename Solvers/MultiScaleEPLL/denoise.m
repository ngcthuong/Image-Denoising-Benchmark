function [cleanI,psnr] = denoise(noiseI,I,models,range,noiseSD,jmp,filters,weights)
% noiseI  - degraded image
% I       - clean image (for computing PSNR)
% models  - cell of GMM models for the different scales
% range   - beta values
% noiseSD - standard deviation of the noise
% jmp     - downsampling factors
% filters - blur filters applied before downsampling
% weights - the weights for the different scales

% params
patchSize = 8;
lambda = patchSize^2/noiseSD^2;

% init clean image
cleanI = noiseI;

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
		upsampled = my_scol2im(cleanZ,patchSize,size(noiseI),jmp(i));

		% backward filter
		comp{i} = imfilter(upsampled,filters{i},'symmetric','conv');
    end

    % conjugate gradient
    b = lambda*noiseI;
    for i=1:length(filters)
        b = b + weights(i)*betas(i)*patchSize^2*comp{i};
    end

    aux_fun = @(x) solveInv_denoise(x,size(noiseI),lambda,weights,betas,patchSize,filters);
    b = b(:);
    [cleanI,~] = cgs(aux_fun,b,1e-15,20);
    cleanI = reshape(cleanI,size(noiseI));

    % psnr
    psnr = [psnr 20*log10(1/std2(cleanI-I))];

%     figure(1); plot(psnr); title('PSNR');
%     drawnow;
end




