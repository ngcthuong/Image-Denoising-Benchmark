function [] = save_image_result(rec, folder_name, img_name, subrate, rec_psnr, rec_ssim, trial)

link = ['ResultImage\' folder_name];
if ~exist(link, 'dir');
    mkdir(link);
end;

final_name = [img_name, '_Sigma', num2str(subrate), '_PSNR', num2str(rec_psnr), 'dB_SSIM', num2str(rec_ssim) ];
if(trial ~= 0)
    final_name = [final_name '_trial' num2str(trial)];
end
final_name = [final_name '.tif'];
if (max(rec(:)) < 10) 
    rec = rec*256; 
end;
imwrite(uint8(rec), strcat([link '\', final_name]));