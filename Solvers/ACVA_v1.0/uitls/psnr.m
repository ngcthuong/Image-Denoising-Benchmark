function res = psnr(hat, star, std)

    if nargin < 3
        std = std2(star);
    end

    res = 10 * ...
          log(std^2 / mean(reshape((hat - star).^2,1,[]))) ...
          / log(10);

end
