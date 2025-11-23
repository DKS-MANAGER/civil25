% Autocorrelation analysis for streamflow (120-day series)
clear; close all; clc;
rng(42);

N = 120; % days
maxLag = 30;

% Parameters (baseline)
alpha = 5; beta = 3; % gamma
prob_rain = 0.4; % p
d = 0.3; c = 0.5; a = 0.6; S_max = 500; S_init = 100;

% Local gamma sampler function must exist: gamrnd_local
if exist('gamrnd_local','file') ~= 2
    error('gamrnd_local.m not found in path.');
end

%% Helper: simulate streamflow implemented in local function 'simulate_streamflow' below

% Generate baseline and compute
rain_baseline = generate_rain(N, prob_rain, alpha, beta);
Q_baseline = simulate_streamflow(rain_baseline, d, c, a, S_max, S_init);

% Part (a): autocorrelation baseline
[acf_baseline, lags] = sample_acf(Q_baseline, maxLag);
threshold = 1.96 / sqrt(N);
sig_lags_baseline = find(abs(acf_baseline(2:end)) > threshold); % exclude lag 0
if isempty(sig_lags_baseline)
    max_sig_lag_baseline = 0;
else
    max_sig_lag_baseline = max(sig_lags_baseline);
end
fprintf('Baseline (p=0.4, c=0.5): max significant lag = %d (95%% CI ±%.3f)\n', max_sig_lag_baseline, threshold);

% Save baseline plot
figure('Position',[100,100,600,400]);
stem(lags, acf_baseline, 'filled'); hold on;
plot([0 maxLag], [threshold threshold], 'r--');
plot([0 maxLag], [-threshold -threshold], 'r--');
xlabel('Lag (days)'); ylabel('Autocorrelation');
title(sprintf('Autocorrelation of streamflow (baseline p=%.2f, c=%.2f)', prob_rain, c));
grid on;
saveas(gcf, 'acf_baseline.png');
close(gcf);

%% Part (b): p = 0.6
prob_rain_b = 0.6;
rain_b = generate_rain(N, prob_rain_b, alpha, beta);
Q_b = simulate_streamflow(rain_b, d, c, a, S_max, S_init);
[acf_b, ~] = sample_acf(Q_b, maxLag);
sig_lags_b = find(abs(acf_b(2:end)) > threshold);
if isempty(sig_lags_b)
    max_sig_lag_b = 0;
else
    max_sig_lag_b = max(sig_lags_b);
end
fprintf('p=0.6 (c=0.5): max significant lag = %d\n', max_sig_lag_b);

figure('Position',[100,100,600,400]);
stem(lags, acf_b, 'filled'); hold on;
plot([0 maxLag], [threshold threshold], 'r--');
plot([0 maxLag], [-threshold -threshold], 'r--');
xlabel('Lag (days)'); ylabel('Autocorrelation');
title(sprintf('Autocorrelation of streamflow (p=%.2f, c=%.2f)', prob_rain_b, c));
grid on;
saveas(gcf, 'acf_p06.png');
close(gcf);

%% Part (c): c = 0.1
c_c = 0.1;
rain_c = generate_rain(N, prob_rain, alpha, beta);
Q_c = simulate_streamflow(rain_c, d, c_c, a, S_max, S_init);
[acf_c, ~] = sample_acf(Q_c, maxLag);
sig_lags_c = find(abs(acf_c(2:end)) > threshold);
if isempty(sig_lags_c)
    max_sig_lag_c = 0;
else
    max_sig_lag_c = max(sig_lags_c);
end
fprintf('c=0.1 (p=0.4): max significant lag = %d\n', max_sig_lag_c);

figure('Position',[100,100,600,400]);
stem(lags, acf_c, 'filled'); hold on;
plot([0 maxLag], [threshold threshold], 'r--');
plot([0 maxLag], [-threshold -threshold], 'r--');
xlabel('Lag (days)'); ylabel('Autocorrelation');
title(sprintf('Autocorrelation of streamflow (p=%.2f, c=%.2f)', prob_rain, c_c));
grid on;
saveas(gcf, 'acf_c01.png');
close(gcf);

fprintf('Plots saved: acf_baseline.png, acf_p06.png, acf_c01.png\n');

%% Local functions
function rain = generate_rain(N, p, alpha, beta)
    rain = zeros(1,N);
    for t=1:N
        if rand() < p
            rain(t) = gamrnd_local(alpha, beta);
        else
            rain(t) = 0;
        end
    end
end

function [streamflow, storage] = simulate_streamflow(rainfall_series, d, c, a, S_max, S_initial)
    n_days = length(rainfall_series);
    storage = zeros(1, n_days + 1);
    surface_runoff = zeros(1, n_days);
    baseflow = zeros(1, n_days);
    infiltration = zeros(1, n_days);
    streamflow = zeros(1, n_days);
    storage(1) = S_initial;
    for t = 1:n_days
        surface_runoff(t) = d * (storage(t) / S_max) * rainfall_series(t);
        baseflow(t) = c * storage(t);
        infiltration(t) = a * rainfall_series(t);
        streamflow(t) = surface_runoff(t) + baseflow(t);
        new_storage = storage(t) + infiltration(t) - baseflow(t);
        storage(t + 1) = max(0, min(new_storage, S_max));
    end
    storage = storage(2:end);
end

function [acf, lags] = sample_acf(x, maxLag)
    x = x(:);
    N = length(x);
    mu = mean(x);
    denom = sum((x - mu).^2);
    acf = zeros(maxLag+1,1);
    for k = 0:maxLag
        if k==0
            acf(k+1) = 1;
        else
            acf(k+1) = sum( (x(1:N-k)-mu) .* (x(k+1:N)-mu) ) / denom;
        end
    end
    lags = (0:maxLag)';
end
