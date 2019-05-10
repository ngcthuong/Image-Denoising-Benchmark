function [psnr, mse] = calPSNR(RGB_rec, RGB_org)

	[H, W, D] = size(RGB_rec);
	mse = sum((RGB_org(:)-RGB_rec(:)).^2)/(H*W*D);
	psnr = 10*log10(255*255/mse);

	% num_pixels = size(RGB_rec, 1) * size(RGB_rec, 2);
    % 
    % % MSE calculation        
    %     err      = RGB_org - RGB_rec;
    %     sse      = sum(sum(err(:).^2)) / 3;
	% 	mse      = sse / num_pixels;              
	% 	
    % % PSNR calculation
    %     psnr     = 10*log10(255^2 ./ mse);
end