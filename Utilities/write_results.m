function write_results(file_name, image, nSig, psnr, ssim, fsim, deltaEInput, trial, inpsnr, inssim, infsim, indeltaEInput)

fid = fopen(file_name, 'a+');
[n, m, c]  = (size(image));

if trial == 0 % take the average
    psnr = mean(psnr);
    ssim = mean(ssim);
    deltaEInput = mean(deltaEInput);
    fsim        = mean(fsim); 
end;


fprintf(fid, '            %6.0f %6.0f %6.0f %6.2f %6.3f %6.3f %6.3f %6.2f %6.3f %6.3f %6.3f' , ... 
        m, trial, nSig, psnr, ssim , fsim, deltaEInput, inpsnr, inssim, infsim, indeltaEInput);
fprintf(fid, '\n');
fclose(fid);
