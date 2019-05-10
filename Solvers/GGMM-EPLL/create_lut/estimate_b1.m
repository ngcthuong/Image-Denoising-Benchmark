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

function b1 = estimate_b1(nu, la)

% Backup warning state and turn off some of them
warnback = warning;
warning('off', 'MATLAB:integral:MinStepSize');
warning('off', 'MATLAB:integral:NonFiniteValue');
warning('off', 'MATLAB:integral:MaxIntervalCountReached');
opts = optimset('display','off');
if length(la) > 1
    f = zeros(size(la));
    for i = 1:size(la, 2)
        b1(:, i) = estimate_b1(nu, la(i));
    end
    return
end

if la == 0
    b1 = + log(sqrt(2 * pi));
    return;
end
alpha = la * exp((gammaln(1./nu) - gammaln(3./nu))/2);

b1 = zeros(size(nu));
s = zeros(size(nu));
for k = 1:length(nu)
    %% -log of the function to be integrated
    g = @(t) t.^2 ./ 2 + exp(nu(k) * (log(abs(t)) - log(alpha(k))));

    %% locate the argmax c and the max s(k)
    c = 0;
    c = fminsearch(g, c); % because gg_shrinkage is approximative
    s(k) = g(c);

    %% find a range that capture enough of the tail
    w = 20*min(1, la);
    t = linspace(c-w, c+w, 1001);
    gt = g(t);
    [s(k), i] = min(gt);
    c = t(i);
    while min(-gt + s(k)) > -14
        w = 2 * w;
        t = linspace(c-w, c+w, 1001);
        gt = g(t);
        [s(k), i] = min(gt);
        c = t(i);
    end

    %% This code was used for debugging
    if 1 == 0
        h = figure;
        subplot(1, 2, 1)
        plot(t, (-g(t) + s(k)))
        xlim([c-w, c+w]);
        subplot(1, 2, 2)
        plot(t, exp(-g(t) + s(k)))
        xlim([c-w, c+w]);
        waitfor(h)
    end

    %% Numerical integration here
    % s(k) is used to ensure that exp(-...) is not
    % zero for all values of x in the range.
    % we are guarantee that it will be 1 at least at c.
    % it will be added back after the loop.
    % We use that adding a constant in exp is a
    % mulitplication.
    q1(k) = integral(@(x) exp(-g(x) + s(k)), ...
                 c-w, c+w, 'AbsTol', 1e-12);
    q2(k) = integral(@(x) x.^2 .* exp(-g(x) + s(k)), ...
                 c-w, c+w, 'AbsTol', 1e-12);
end
b1 = -log(2) + log(1-q2./q1);

% Restitute initial warning state
warning(warnback);
