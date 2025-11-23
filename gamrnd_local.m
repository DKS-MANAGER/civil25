function r = gamrnd_local(k, theta, varargin)
% GAMRND_LOCAL  Local gamma random number generator (shape k, scale theta)
% Usage:
%   r = gamrnd_local(k, theta)
%   r = gamrnd_local(k, theta, m, n)  -> returns m-by-n matrix
%
% Uses Marsaglia & Tsang method and supports scalar or matrix outputs.

% Determine output size
if isempty(varargin)
    outSize = [1 1];
else
    outSize = cell2mat(varargin);
end

num = prod(outSize);
r = zeros(outSize);

for idx = 1:num
    kk = k;
    if kk <= 0
        error('Shape parameter k must be > 0');
    end
    if kk == 1
        sample = -log(rand());
    elseif kk < 1
        % Use boost: Gamma(k) = Gamma(k+1) * U^(1/k)
        kk1 = kk + 1;
        d = kk1 - 1/3;
        c = 1 / sqrt(9*d);
        while true
            x = randn();
            v = (1 + c*x)^3;
            if v > 0
                u = rand();
                if u < 1 - 0.0331*(x^4) || log(u) < 0.5*x^2 + d*(1 - v + log(v))
                    g = d * v;
                    break;
                end
            end
        end
        u2 = rand();
        sample = g * (u2^(1/kk));
    else
        d = kk - 1/3;
        c = 1 / sqrt(9*d);
        while true
            x = randn();
            v = (1 + c*x)^3;
            if v > 0
                u = rand();
                if u < 1 - 0.0331*(x^4) || log(u) < 0.5*x^2 + d*(1 - v + log(v))
                    sample = d * v;
                    break;
                end
            end
        end
    end
    r(idx) = theta * sample;
end
end