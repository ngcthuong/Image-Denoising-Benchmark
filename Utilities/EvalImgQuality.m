function [out] = EvalImgQuality2(ImgRec, ImgOrg, Qtype, varargin)

if ~isempty(varargin)
    dynamic = varargin{1};
else
    dynamic = 1;
end
if strcmp(Qtype, 'PSNR')
    MSE = sum(sum(sum((ImgRec - ImgOrg).^2)));
    MSE= MSE/(numel(ImgRec));
    out = 10 * log10(dynamic^2/MSE);
    
elseif strcmp(Qtype, 'Delta2000')
    
    %wp = whitepoint('d65');
    C = makecform('srgb2lab');
    ImgRec = ImgRec/dynamic;
    ImgOrg = ImgOrg/dynamic;
    [y, x, z] = size(ImgRec);
    
    Ild = applycform(ImgRec,C);
    Jld = applycform(ImgOrg,C);
    
    Ilds = [reshape(Ild(:,:,1),[y*x 1]) reshape(Ild(:,:,2),[y*x 1])  reshape(Ild(:,:,3),[y*x 1]) ];
    Jlds = [reshape(Jld(:,:,1),[y*x 1]) reshape(Jld(:,:,2),[y*x 1])  reshape(Jld(:,:,3),[y*x 1]) ];
    
    d2k = deltaE2000(Ilds, Jlds);
    out= mean(d2k(:));    
end