function write_results2(file_name, image, nSig, psnr1, ssim1, psnr2, ssim2, trial, inpsnr, inssim, psnr3, ssim3)

fid = fopen(file_name, 'a+');
[n, m, c]  = (size(image));

if trial == 0 % take the average
    psnr1 = mean(psnr1);
    ssim1 = mean(ssim1);
    ssim2 = mean(ssim2);
    psnr2        = mean(psnr2); 
end;


fprintf(fid, '            %6.0f %6.0f %6.0f %6.2f %6.3f %6.3f %6.3f %6.2f %6.3f %6.3f %6.3f' , ... 
        m, trial, nSig, psnr1, ssim1 , psnr2, ssim2, inpsnr, inssim, psnr3, ssim3);
fprintf(fid, '\n');
fclose(fid);
