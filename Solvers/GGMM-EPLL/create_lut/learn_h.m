% % Undocumented.
%
% C-A. Deledalle, S. Parameswaran, and T. Q. Nguyen, "Image
% restoration with generalized Gaussian mixture model patch
% priors", arXiv.
%
% If you use any part of this software package please cite the above
% publication
%
%
% License details as in license.txt
% ________________________________________

clear all
close all

cd(fileparts(mfilename('fullpath')));
addpath('..');
addpathrec('..');

L = 100;
N = 100;

la = logspace(log10(0.001), log10(1000), L);
nu = logspace(log10(0.3), log10(2), N);
h = zeros(L, N);
for k = 1:L
    fprintf('%.2f%%', (k-1)/(L-1));
    for l = 1:length(nu)
        a1 = gg_a1(la(k), nu(l));
        b1 = gg_b1(la(k), nu(l));
        a2 = gg_a2(la(k), nu(l));
        b2 = gg_b2(la(k), nu(l));

        gamma = gg_gamma(la(k), nu(l), 'direct');
        if a1 == a2
            c = 2;
        else
            c = (b2-b1)/(a1-a2);
        end
        x = logspace(-3, 2*log10(exp(c))+3, 101);
        y = log(gg_discrepancy(x, 1, la(k), nu(l), 'mode', 'direct') - gamma);
        lx = log(x);
        f1 = a1 .* lx + b1;
        f2 = a2 .* lx + b2;

        y = real(y);
        for i = 100:-1:1
            if y(i) > y(i+1)
                y(i) = f1(i);
            end
            if y(i) > y(i+1)
                y(i) = y(i+1);
            end
        end
        for i = 2:101
            if isinf(y(i))
                y(i) = f2(i);
            end
        end

        if nu(l) <= 2
            f = @(h) bsxfun(@minus, f2', ...
                            bsxfun(@times, h, ...
                                   log(1+exp((f2-f1)' * (1./h)))));
        else
            f = @(h) bsxfun(@plus, f1', ...
                            bsxfun(@times, h, ...
                                   log(1+exp((f2-f1)' * (1./h)))));
        end
        cost = @(h) sum(bsxfun(@minus, f(h), y').^2, 1);
        h(k, l) = exp(fminsearch(@(h) cost(exp(h)), log(0.1)));

        continue
        fancyfigure
        plot(lx, y, 'Linewidth', 2);
        hold on
        plot(lx, f2, '--');
        plot(lx, f1, '--');
        plot(lx, f(h(k, l)), '--');
        waitfor(gcf)
    end
end

hn = h;
l = length(nu);
for k = 1:length(la)
    hn(k,l) = h(k,l-1);
end
hn(1:8, 99) = hn(1:8, 98);
hn(1:8, 100) = hn(1:8, 98);

close all
fancyfigure
imagesc(log10(nu), log10(la), h)
set(gca, 'XTickMode', 'manual');
set(gca, 'YTickMode', 'manual');
xlabel('$\log_{10} \lambda$');
ylabel('$\log_{10} \nu$');
title('$a_1$');
axis square

return

log_la = log(la);
log_nu = log(nu);
save('h_lookup_table_100x100.mat', 'log_la', 'log_nu', 'h');
