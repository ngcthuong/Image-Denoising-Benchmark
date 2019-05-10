clear all;  close all;  path(path,genpath(pwd));
imgSize     = 256;              % 
testImg     = [1];           %  
testSig     = [ 30 50 70];    % 
recMode     = { 'BM3D'      'WNNM'      'GSRC',     'AST-NLS',  'PGPD', ...
                'MSEPLL',   'DnCNN',    'SSC_GSM'   'ACPT'   ,  'TWSC'  ...
                'NCSR'};
recMode_id  = 11;

inPar.isShow= 0;
nbrTrial    = 1;
strNote     = ['_' recMode{recMode_id} '_20190510'];
disp( ['Start Simulation '  'Denoising - ' strNote] );

for imId = 1:1:length(testImg)
    % read out the image
    [ImgOrg, inPar.imgName] = testImage(imgSize, testImg(imId));
    ImgOrg        	= double(ImgOrg);
    [n1, n2, n3]  	= size(ImgOrg);
    
    % Define the output file 
    saveFolderText    = ['ResultText' num2str(n1) '\' ];   
    if ~exist(saveFolderText, 'dir');  mkdir(saveFolderText);   end;    
    fileNameSave  	  = [saveFolderText inPar.imgName  strNote ];
    write_info([fileNameSave  '.txt'], ['%' inPar.imgName ]);  
    write_info([fileNameSave  '.txt'], ['%             size  trial    Sub  PSNR   SSIM   FSIM   Delt  nPSNR   nSSIM  nFSIM  nDelta']);
    fileNameSaveAll   = [saveFolderText 'all_' strNote ];
    
    for sigId = 1:1:length(testSig)
        inPar.nSig = testSig(sigId);  display(['   ' num2str(imgSize) '_' inPar.imgName '_' strNote  '_sigma' num2str(inPar.nSig)] );        
        for trial = 1:1:nbrTrial
            % noisy image            
            [ImgNoise]  = addNoise(ImgOrg, inPar.nSig, trial);
             % denoise image
            [ImgRec, outPar]= Denoiser(ImgNoise, ImgOrg, recMode{recMode_id}, inPar);            
            % result plot
            [inPSNR(trial), inSSIM(trial)] = EvalImgQuality2(ImgNoise, ImgOrg);
            [outPSNR(trial), outSSIM(trial)] = EvalImgQuality2(ImgRec, ImgOrg);
            disp(['      trial: ' num2str(trial) ', Noisy = ', num2str(inPSNR(trial))  'dB, Denoised = ' num2str(outPSNR(trial)) 'dB']);
        end
        % write average results        
        write_results([fileNameSave  '.txt'], ImgOrg, inPar.nSig, mean(outPSNR), mean(outSSIM), mean(0), ...
                            mean(0), trial, mean(inPSNR), mean(inSSIM), mean(0), mean(0))
        save_image_result(ImgRec, [recMode{recMode_id}],  inPar.imgName, inPar.nSig, mean(outPSNR), mean(outSSIM), 0);   
        % save all resutls to one file
        save([fileNameSaveAll  '_nSig' num2str(inPar.nSig) '.mat'], 'ImgOrg', 'inPSNR', 'inSSIM', ...
                                        'outPSNR', 'outSSIM', 'inPar', 'outPar', 'ImgRec');        
    end % ending sigma
end % ending test image
disp('SIMULATION END!!!');