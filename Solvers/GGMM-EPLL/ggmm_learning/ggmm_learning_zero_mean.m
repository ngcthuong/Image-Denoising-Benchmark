function [ gmm, energy, gamma ] = ggmm_learning_zero_mean(x, K, R, gmm, nu0)
% % Function Name: ggmm_learning_zero_mean
%
%
% Undocummented.

% Citation:
% If you use this code please cite:
%
% C-A. Deledalle, S. Parameswaran, and T. Q. Nguyen, "Image
% restoration with generalized Gaussian mixture model patch
% priors", arXiv.
%
% License details as in license.txt
% ________________________________________

eps = 10^-6;
[n, Q] = size(x);
if ~exist('R', 'var')
    R = 10;
end

% Initialization
fprintf('Start learning: K=%d N=%d\n', K, Q);
if ~exist('gmm', 'var') || isempty(gmm)
   error('Initialization required');
else
    nu          = gmm.nu;
    Sigma       = gmm.Sigma;
    SSigma      = gmm.SSigma;
    USigma      = gmm.USigma;
    iSigma      = gmm.iSigma;
    cSigma      = gmm.cSigma;
    logdetSigma = gmm.logdetSigma;
    weight      = gmm.weight;
end
if exist('nu0', 'var')
  for k =1:K
    nu{k}(:) = nu0;
  end
end

% Lookuptable
if exist('Finv_0_075_1024.mat', 'file')
    load('Finv_0_075_1024.mat', 'Finv');
else
    K = 1024;
    x = linspace(0, 200, K);
    F = @(x) exp(2 * gammaln(2 ./ x) - gammaln(3 ./ x) - gammaln(1 ./ x));
    Fx = F(x);

    y_list = linspace(0, 0.75, K);
    Finv = zeros(K, 1);
    for i = 1:K
        %i
        y = y_list(i);
        xmin = 0;
        xmax = 200;
        while abs((xmax - xmin)) > 1e-10
            xc = (xmax + xmin) / 2;
            if F(xc) > y
                xmax = xc;
            else
                xmin = xc;
            end
        end
        Finv(i) = xc;
    end

    subplot(1,2,1);
    plot(x, Fx);
    subplot(1,2,2);
    plot(y_list, Finv);
    save('Finv_0_075_1024.mat', 'Finv');
end

% eta function
eta = @(nu) exp((gammaln(3./nu)-gammaln(1./nu))/2);

gamma = zeros(Q, K);
llk   = zeros(Q, K);
energy = zeros(R, 1);
r = 1;
while r <= R
    % E-Step
    redo = 1;
    while redo
        redo = 0;
        for k = 1:K
            a         = sqrt(SSigma{k}) ./ eta(nu{k});
            diff      = USigma{k}' * x;
            diff      = bsxfun(@rdivide, diff, a);
            llk(:, k) = ...
                sum(bsxfun(@power, abs(diff), nu{k}), 1) + ...
                - sum(log(nu{k})) + ...
                + sum(log(2 * a) + gammaln(1 ./ nu{k}));
        end
        nllk  = bsxfun(@minus, llk, min(llk, [], 2));
        gamma = bsxfun(@times, exp(-nllk), weight);
        gamma = bsxfun(@times, gamma, 1./sum(gamma, 2));

        fprintf('Iteration %d:', r);
        fprintf(' sc %.2f bc %.2f', min(sum(gamma)), max(sum(gamma)));

        if r > 1
            t = 2.1;
            % Resample unused Gaussian k by splitting another one
            for k = find(sum(gamma) <= t)
                redo = 1;
                % Choose Gaussian i to split with proba sum(gamma)
                cs   = cumsum(sum(gamma));
                cs   = cs / cs(end);
                csi  = [0 cs(1:(end-1))];
                i    = k;
                while sum(gamma(:,i)) <= t
                    p = rand;
                    i = find(csi <= p & p < cs);
                    i = i(1);
                end
                dir            = cSigma{i} * randn(n, 1);
                nu{k}          = nu{i};
                nu{i}          = nu{i};
                Sigma{k}       = Sigma{i};
                SSigma{k}      = SSigma{i};
                USigma{k}      = USigma{i};
                iSigma{k}      = iSigma{i};
                cSigma{k}      = cSigma{i};
                logdetSigma{k} = logdetSigma{i};
                weight(k)      = weight(i) / 2;
                weight(i)      = weight(i) / 2;
            end
            if redo
                fprintf(' (Resample)\n');
            end
        end
    end

    % M-Step
    Z = sum(gamma);
    [~, idxw] = sort(weight, 'descend');
    for kp = 1:K
        k = idxw(kp);
        if Z(k) == 0
            continue;
        end
        weight(k)           = Z(k) / sum(Z);
        xg                  = bsxfun(@times, x, gamma(:, k)');
        Sigma{k}            = xg * x' / Z(k);
        [E, D]              = mysvd(Sigma{k});
        if D(1) == 0
            D(1) = 1;
        end
	D(D == 0) = min(D(D > 0));

        Exag                = bsxfun(@times, abs(E' * x), gamma(:, k)');
        D1                  = sum(Exag, 2) / Z(k);
        Fnu                 = D1.^2 ./ D;
        Fnu(Fnu > 0.75)     = 0.75;
        idx                 = 1 + round(Fnu / 0.75 * 1023);
        nu{k}               = Finv(idx);
        nu{k}(nu{k} > 2)    = 2;
        nu{k}(nu{k} < 0.3)  = 0.3;
	if exist('nu0', 'var')
	  nu{k}(:) = nu0;
	end

        D(D < eps*D(1))     = eps*D(1);
        SSigma{k}           = D;
        USigma{k}           = E;
        iSigma{k}           = E * diag(1./D)       * E';
        Sigma{k}            = E * diag(D)          * E';
        cSigma{k}           = E * diag(sqrt(D))    * E';
        logdetSigma{k}      = sum(log(D));
    end

    % Energy
    energy(r) = ...
        mean(sum(bsxfun(@times, gamma, log(weight)))) + ...
        mean(sum(gamma .* (-llk)));

    fprintf(' e %.2f H %.2f', energy(r), ...
            -sum(weight.*log(weight))/log(K));
    fprintf('\n');

    r = r + 1;
end
clear llk

% Sort Gaussians by weights
[~, idx]          = sort(weight, 'descend');
gmm.nu            = nu(idx);
gmm.Sigma         = Sigma(idx);
gmm.iSigma        = iSigma(idx);
gmm.cSigma        = cSigma(idx);
gmm.logdetSigma   = logdetSigma(idx);
for k = 1:K
    [U, S]        = mysvd(gmm.Sigma{k});
    gmm.USigma{k} = U;
    gmm.SSigma{k} = S;
end
gmm.weight        = weight(idx);
gmm.logweight     = log(gmm.weight);



function varargout = mysvd(A)
% same as SVD but D is a vector

    [U, D] = svd(A);
    D = diag(D);
    if nargout == 1
        varargout{1} = D;
    end
    if nargout == 2
        varargout{2} = D;
        varargout{1} = U;
    end
