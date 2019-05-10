function Par = ParSetDnCNN(nSig, input, useGPU)

Par.input   = input; 
folderModel = 'model';
noiseSigma  = nSig;  %%% image noise level
showResult  = 1;
Par.useGPU  = useGPU;
pauseTime   = 0;

%%% load [specific] Gaussian denoising model

modelSigma  = min(75,max(10,round(noiseSigma/5)*5)); %%% model noise level
load(fullfile(folderModel,'specifics',['sigma=',num2str(modelSigma,'%02d'),'.mat']));

Par.net = net; 

%%% move to gpu
if useGPU
    Par.net = vl_simplenn_move(net, 'gpu') ;
end
    
%%% convert to GPU
if useGPU
	Par.input = gpuArray(input);
end




