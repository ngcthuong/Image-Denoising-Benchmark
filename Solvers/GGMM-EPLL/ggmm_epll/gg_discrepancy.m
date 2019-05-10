function f = gg_discrepancy(x, si, la, nu, varargin)
% % Function Name: gg_discrepancy
%
%
% Inputs:
%   x           : real value or array
%   si          : Gaussian noise std
%   la          : GGD standard deviation
%   nu          : Shape parameter of the GGD
%   varargin    : refer to retrieve arguments for a list
%
% Outputs:
%   f           : value of the discrepancy function

% Citation:
% If you use this code please cite:
%
% C-A. Deledalle, S. Parameswaran, and T. Q. Nguyen, "Image
% restoration with generalized Gaussian mixture model patch
% priors", arXiv.
%
% License details as in license.txt
% ________________________________________


% Retrieve arguments
options = makeoptions(varargin{:});
mode    = getoptions(options, 'mode', 'lut');

% Reduction
tau = si;
x   = x ./ tau;
si  = si ./ tau;
la  = la ./ tau;

for dummy = 1 % Used to break/quit the switch
    switch mode
        case 'stochastic' % Very poor performance
            alpha = la * exp((gammaln(1/nu) - gammaln(3/nu))/2);
            for k = 1:length(x)
                g = @(t) -abs(t/alpha).^nu;
                w = x(k) + randn(10000, 1);
                g = g(w);
                m = max(g);
                f(k) = log(mean(exp(g - m))) + m;
            end
            f = log(2 * alpha / nu) ...
                + gammaln(1/nu) ...
                - f;
        case 'direct' % Very slow
            if length(x) == 1 && length(la) > 1
                f = zeros(size(la));
                for i = 1:size(la, 1)
                    for j = 1:size(la, 2)
                        f(i, j) = gg_discrepancy(x, 1, la(i, j), nu, ...
                                                 'mode', mode);
                    end
                end
                break;
            end
            if size(x, 1) > 1 && size(x, 2) > 1
                f = zeros(size(x));
                for k = 1:size(x, 2)
                    f(:, k) = gg_discrepancy(x(:, k), 1, la, nu, ...
                                             'mode', mode);
                end
                break;
            end
            % Backup warning state and turn off some of them
            warnback = warning;
            warning('off', 'MATLAB:integral:MinStepSize');
            warning('off', 'MATLAB:integral:NonFiniteValue');
            warning('off', 'MATLAB:integral:MaxIntervalCountReached');
            opts = optimset('display','off');

            if length(la) == 1
                la = repmat(la, size(x));
            end
            if length(nu) == 1
                nu = repmat(nu, size(x));
            end
            if norm(size(la) - size(x)) ~= 0 || norm(size(nu) - size(x)) ~= 0
                error(sprintf(['Dimension issue:\n' ...
                               'x is %d %d\n' ...
                               'la is %d %d\n' ...
                               'nu is %d %d'], ...
                              size(x, 1), size(x, 2), ...
                              size(la, 1), size(la, 2), ...
                              size(nu, 1), size(nu, 2)));
            end
            alpha = la .* exp((gammaln(1./nu) - gammaln(3./nu))/2);
            f = zeros(size(x));
            s = zeros(size(x));
            for k = 1:length(x)
                if la(k) == 0
                    continue;
                end

                %% -log of the function to be integrated
                g = @(t) (x(k)-t).^2 ./ 2 + exp(nu(k) * (log(abs(t)) - log(alpha(k))));

                %% locate the argmax c and the max s(k)
                c = gg_shrinkage(x(k), 1, la(k), nu(k));
                c = fminsearch(g, c); % because gg_shrinkage is approximative
                s(k) = g(c);

                %% find a range that capture enough of the tail
                w = 20*min(1, la(k));
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
                q = integral(@(x) exp(-g(x) + s(k)), ...
                             c-w, c+w, 'AbsTol', 1e-12);
                f(k) = log(q);
                if isnan(f(k)) || isinf(f(k))
                    % the function is probably a dirac of mass 1
                    f(k) = 0;
                end
            end
            f = + log(sqrt(2 * pi)) + ...
                + log(2 * alpha ./ nu) ...
                + gammaln(1./nu) ...
                + s ...
                - f;
            f(la == 0) = ...
                  + log(sqrt(2 * pi)) + ...
                  + x(la == 0).^2 ./ 2;

            % Restitute initial warning state
            warning(warnback);
            break;
        case 'lut' % Fast and good performance
            lut = getoptions(options, 'lut', []);
            if isempty(lut)
                load('params_lookup_table.mat', 'lut')
            end
            [ga, a1, a2, b1, b2, h] = gg_discrepancy_params(la, nu, lut);
            switch 'softplusc'
                case 'relu'
                    lx  = log(abs(x));
                    f1  = a1 .* lx + b1;
                    f2  = a2 .* lx + b2;
                    idx = (nu <= 2);
                    f   = zeros(size(f1));
                    f(idx,:)  = min(f1(idx,:),  f2(idx,:));
                    f(~idx,:) = max(f1(~idx,:), f2(~idx,:));
                    f   = exp(f) + ga;
                case 'softplus'
                    lx  = log(abs(x));
                    f1  = a1 .* lx + b1;
                    f2  = a2 .* lx + b2;
                    di  = (f2 - f1)  ./ h;
                    le  = h .* log(1 + exp(di));
                    f   = (nu <  2) .* (f2 - le) + ...
                          (nu >  2) .* (f1 + le) + ...
                          (nu == 2) .* f1;
                    f(isinf(lx)) = -Inf;
                    f   = exp(f) + ga;
                case 'softplusc'
                    % Same as softplus but written in MEX-C
                    f = ga + exp(gg_discrepancy_core(log(abs(x)), nu, a1, a2, b1, b2, h));
            end
    end
end
f = f + log(tau);