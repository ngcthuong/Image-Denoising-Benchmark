function [Par] = ParSetMSEPLL(nSig)
		% models
        nSig = nSig/255; 
		load GSModel_8x8_200_2M_noDC_zeromean;
		load GMM_high;
		Par.models = {GS,GS,GS};
		% jump size
		Par.jmp = [1,2,4];

		% filters
		G = fspecial('gaussian',11,0.8);
		GG = fspecial('gaussian',11,1.5);
		H = zeros(size(G));
		H((1+size(H,1))/2,(1+size(H,2))/2) = 1;
		H = H - fspecial('gaussian',11,0.6);

		Par.filters = {1,G,GG};

		% weights
		Par.weights = [1,0.05,0.05];
		Par.betas 	= (1/nSig^2) *[1 4 8 16 32 64];
		
