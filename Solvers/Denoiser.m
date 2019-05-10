function [ImgRec, res]= Denoiser(ImgNoise, ImgOrg, recMode, inPar)
res         = [];
% par         = ParSetDenoise(recMode);
nSig    = inPar.nSig; 
isShow  = inPar.isShow; 

switch recMode
   case 'BM3D'
        [res.PSNR, ImgRec] = BM3D(ImgOrg, ImgNoise, nSig );
        ImgRec = ImgRec * 255; 
    case 'WNNM'
        par     = ParSetWNNM(nSig);   
        par.step= 1; 
        ImgRec  = WNNM_DeNoising( ImgNoise, ImgOrg, par );         
    case 'GSRC'
        par     = ParSetGSRC (nSig, ImgOrg);
        par.nim = ImgNoise; 
        par.I   = ImgOrg; 
        ImgRec  = GSRC_Denoising( par, par.Thr);		
	case 'AST-NLS'
		ImgRec  = ast_nls( ImgNoise, nSig);	
	case 'MSEPLL'
		% models	
		par 	= ParSetMSEPLL(nSig);
		[ImgRec, res.PSNR] = denoise(ImgNoise/255, ImgOrg/255, par.models, par.betas, nSig, par.jmp, par.filters, par.weights);
		ImgRec 	= ImgRec * 255; 		
	case 'PGPD'
		[par, model]  	=  ParSetPGPD( nSig );
		par.I 			= ImgOrg/255; 
		par.nim 		= ImgNoise/255; 
		[ImgRec, par]  	=  PGPD_Denoising(par, model);		
	case 'SSC_GSM'
		if nSig < 50  	K = 3; 
        elseif nSig < 100 K = 4; 
		else 			K = 5;   end;
		par 	= ParSetSSC( nSig, K);
		par.nim = ImgNoise; 
        par.I   = ImgOrg; 
		[ImgRec res.PSNR res.SSIM]   =    SSC_GSM_Denoising( par );    

	%% Deep learning based method
	case 'DnCNN'
		useGPU  = 0; 
		par 	= ParSetDnCNN(nSig, ImgNoise/255, useGPU);
		res 	= simplenn_matlab(par.net, par.input);
		
		if par.useGPU
			ImgRec = gather(ImgRec);
			Par.input  = gather(Par.input);
        end
        ImgRec 	= par.input - res(end).x;
        ImgRec  = ImgRec * 255; 
     %% Your Method
    case 'oAR_GPA_GSRC'
        par     = ParSetGPA_GSRC (nSig, ImgOrg);
        par.nim = ImgNoise; 
        par.I   = ImgOrg; 
        ImgRec  = oAR_GPA_GSRC_Denoising( par, par.Thr );
    case 'oTopN_GPA_GSRC'
        par     = ParSetGPA_GSRC (nSig, ImgOrg);
        par.nim = ImgNoise; 
        par.I   = ImgOrg; 
        ImgRec  = oTopN_GPA_GSRC_Denoising( par, par.Thr );
    case 'oAR_GPA_WNNM'
        RefImg   = ImgOrg; 
        par2     = ParSetWNNM(nSig);   
        par2.step= 1;         
        [ImgRec, par2]  = sAR_GPA_WNNM( ImgNoise, ImgOrg, RefImg, par2 ); 
    case 'oTopN_GPA_WNNM'
        RefImg   = ImgOrg; 
        par      = ParSetWNNM(nSig);   
        par.step = 1;         
        [ImgRec, par2]  = topN_GPA_WNNM( N_Img, O_Img, R_Img, par )
    case 'sAR_GPA_WNNM_DnCNN'
        % Reference with 
        useGPU  = 0; 
		par 	= ParSetDnCNN(nSig, ImgNoise/255, useGPU);
		res2 	= simplenn_matlab(par.net, par.input);		
		if par.useGPU
			ImgRec = gather(ImgRec);
			Par.input  = gather(Par.input);
        end
        ImgRec 	= par.input - res2(end).x;
        RefImg  = ImgRec * 255; 
        disp([inPar.imgName '_Ref_DnCNN_' num2str(calPSNR(RefImg, ImgOrg)) 'dB']);
        % GPA - WNNM 
        par2     = ParSetWNNM(nSig);   
        par2.step= 1;         
        [ImgRec, par2]  = sAR_GPA_WNNM( ImgNoise, ImgOrg, RefImg, par2 );  
        
     case 'sAR_GPA_WNNM_DnCNN_Fuse'
        % Reference with 
        useGPU  = 0; 
		par 	= ParSetDnCNN(nSig, ImgNoise/255, useGPU);
		res2 	= simplenn_matlab(par.net, par.input);		
		if par.useGPU
			ImgRec = gather(ImgRec);
			Par.input  = gather(Par.input);
        end
        ImgRec      = par.input - res2(end).x;
        RefImg      = ImgRec * 255; 
        disp([inPar.imgName '_Ref_DnCNN_' num2str(calPSNR(RefImg, ImgOrg)) 'dB']);
        res.PSNR    = calPSNR(RefImg, ImgOrg); 
        res.RefImg  =RefImg; 
        % GPA - WNNM 
        par2        = ParSetWNNM(nSig);   
        par2.step   = 1;       
        par2.ctqSize= 3;
        par2.thresh  = 0.12; 
        [ImgRec, par2]  = sAR_GPA_WNNM_Fuse( ImgNoise, ImgOrg, RefImg, par2 ); 
        res.par = par; 
     case 'sAR_GPA_WNNM_DnCNN_Fuse2'
        % Reference with 
        useGPU  = 0; 
		par 	= ParSetDnCNN(nSig, ImgNoise/255, useGPU);
		res2 	= simplenn_matlab(par.net, par.input);		
		if par.useGPU
			ImgRec = gather(ImgRec);
			Par.input  = gather(Par.input);
        end
        ImgRec      = par.input - res2(end).x;
        RefImg      = ImgRec * 255; 
        disp([inPar.imgName '_Ref_DnCNN_' num2str(calPSNR(RefImg, ImgOrg)) 'dB']);
        res.PSNR    = calPSNR(RefImg, ImgOrg); 
        res.RefImg  =RefImg; 
        % GPA - WNNM 
        par2        = ParSetWNNM(nSig);   
        par2.step   = 1;       
        par2.ctqSize= 6;
        par2.win    = 16; 
        par2.thresh = 0.12; 
        [ImgRec, par2]  = sAR_GPA_WNNM_Fuse2( ImgNoise, ImgOrg, RefImg, par2 ); 
        res.par = par; 
     case 'sAR_GPA_WNNM_DnCNN_Fuse3'
        % Reference with 
        useGPU  = 0; 
		par 	= ParSetDnCNN(nSig, ImgNoise/255, useGPU);
		res2 	= simplenn_matlab(par.net, par.input);		
		if par.useGPU
			ImgRec = gather(ImgRec);
			Par.input  = gather(Par.input);
        end 
        ImgRec      = par.input - res2(end).x;
        RefImg      = ImgRec * 255; 
        disp([inPar.imgName '_Ref_DnCNN_' num2str(calPSNR(RefImg, ImgOrg)) 'dB']);
        res.PSNR    = calPSNR(RefImg, ImgOrg); 
        res.RefImg  =RefImg; 
        % GPA - WNNM 
        par2        = ParSetWNNM(nSig);   
        par2.step   = 1;  
        par2.thresh = 0.12; 
        [ImgRec, par2]  = sAR_GPA_WNNM_Fuse3( ImgNoise, ImgOrg, RefImg, par2 ); 
        res.par = par; 
    case 'topN_GPA_WNNM_DnCNN'
        useGPU  = 0; 
		par 	= ParSetDnCNN(nSig, ImgNoise/255, useGPU);
		res 	= simplenn_matlab(par.net, par.input);		
		if par.useGPU
			ImgRec = gather(ImgRec);
			Par.input  = gather(Par.input);
        end
        ImgRec 	= par.input - res(end).x;
        RefImg  = ImgRec * 255; 
        disp([inPar.imgName '_Ref_DnCNN_' num2str(csnr(RefImg, ImgOrg, 0, 0)) 'dB']);
        % GPA - WNNM 
        par2     = ParSetWNNM(nSig);   
        par2.step= 1;         
        par2.noPatch = 30; 
        [ImgRec, par2]  = topN_GPA_WNNM( ImgNoise, ImgOrg, RefImg, par2 );  
end;