function [ImgNoise] = addNoise(ImgOrg, nSig, trial)
% Functio [ImgNoise] = addNoise(ImgOrg, sigma, trial)
% Function to add noise to clean image 
%   + Input
%       - 

[n, m, c] = size(ImgOrg); 
noise_filename = ['Noisy/Gauss' num2str(c) '_nSig' num2str(nSig) '_trial' num2str(trial)];
if(exist(noise_filename, 'file'));
    load(noise_filename);
else    
    Noise = zeros(n, m, c); 
    for i = 1:1:c
        Noise(:, :, i) = randn(n, m);
    end;
    save(noise_filename, 'Noise');
end;

% noisy image
ImgNoise      = ImgOrg + nSig.* Noise;
