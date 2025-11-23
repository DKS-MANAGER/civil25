% Stochastic Streamflow Modeling - MATLAB Implementation
% This code implements a complete solution for the monsoon streamflow modeling problem

%% Clear workspace and set random seed
clear all; close all; clc;
rng(42); % Set random seed for reproducibility

%% Parameters
n_simulations = 1000;
n_days = 120;
prob_rain = 0.4;
alpha = 5;
beta = 3;
d = 0.3;
c = 0.5;
a = 0.6;
S_max = 500;
S_initial = 100;

%% Part (a): Generate rainfall time-series
fprintf('Part (a): Generating 1000 rainfall time-series of 120 days each\n');

rainfall_data = zeros(n_simulations, n_days);

for sim = 1:n_simulations
    for day = 1:n_days
        % Bernoulli trial for rain occurrence
        if rand() < prob_rain
            % Generate rainfall from Gamma distribution
            rainfall_data(sim, day) = gamrnd(alpha, beta);
        else
            rainfall_data(sim, day) = 0;
        end
    end
end

% Verify statistical properties
rainy_days = rainfall_data > 0;
observed_prob_rain = mean(rainy_days(:));
rainy_values = rainfall_data(rainfall_data > 0);

fprintf('Observed probability of rain: %.3f (Expected: 0.4)\n', observed_prob_rain);
fprintf('Mean rainfall on rainy days: %.2f mm\n', mean(rainy_values));
fprintf('Std rainfall on rainy days: %.2f mm\n', std(rainy_values));
fprintf('Theoretical mean: %.2f mm\n', alpha * beta);
fprintf('Theoretical std: %.2f mm\n', sqrt(alpha) * beta);

%% Part (b): Streamflow modeling function
function [streamflow, storage, surface_runoff, baseflow, infiltration] = ...
    simulate_streamflow(rainfall_series, d, c, a, S_max, S_initial)
    
    n_days = length(rainfall_series);
    
    % Initialize arrays
    storage = zeros(1, n_days + 1);
    surface_runoff = zeros(1, n_days);
    baseflow = zeros(1, n_days);
    infiltration = zeros(1, n_days);
    streamflow = zeros(1, n_days);
    
    % Set initial storage
    storage(1) = S_initial;
    
    for t = 1:n_days
        % Calculate surface runoff: sr_t = d * (S_{t-1} / S_max) * x_t
        surface_runoff(t) = d * (storage(t) / S_max) * rainfall_series(t);
        
        % Calculate baseflow: bf_t = c * S_{t-1}
        baseflow(t) = c * storage(t);
        
        % Calculate infiltration: inf_t = a * x_t
        infiltration(t) = a * rainfall_series(t);
        
        % Calculate streamflow: Q_t = sr_t + bf_t
        streamflow(t) = surface_runoff(t) + baseflow(t);
        
        % Update storage: S_t = S_{t-1} + inf_t - bf_t
        new_storage = storage(t) + infiltration(t) - baseflow(t);
        storage(t + 1) = max(0, min(new_storage, S_max));
    end
    
    % Return only the storage values for days 1 to n_days
    storage = storage(2:end);
end

% Test the streamflow model
[test_streamflow, test_storage, ~, ~, ~] = ...
    simulate_streamflow(rainfall_data(1,:), d, c, a, S_max, S_initial);

fprintf('\nPart (b): Streamflow modeling completed\n');
fprintf('Maximum streamflow in test simulation: %.2f mm\n', max(test_streamflow));
fprintf('Mean daily streamflow in test simulation: %.2f mm\n', mean(test_streamflow));

%% Part (c): Calculate maximum streamflow for all simulations
fprintf('\nPart (c): Running all 1000 simulations...\n');

max_streamflows_beta3 = zeros(1, n_simulations);

for sim = 1:n_simulations
    [streamflow, ~, ~, ~, ~] = simulate_streamflow(rainfall_data(sim,:), d, c, a, S_max, S_initial);
    max_streamflows_beta3(sim) = max(streamflow);
    
    if mod(sim, 200) == 0
        fprintf('Completed %d simulations\n', sim);
    end
end

fprintf('Mean maximum streamflow (β=3): %.2f mm\n', mean(max_streamflows_beta3));
fprintf('Std maximum streamflow (β=3): %.2f mm\n', std(max_streamflows_beta3));

% Plot probability distribution for β=3
figure('Position', [100, 100, 1200, 400]);

subplot(1,3,1);
histogram(max_streamflows_beta3, 30, 'Normalization', 'probability', 'FaceColor', 'blue', 'FaceAlpha', 0.7);
title(sprintf('β=3, c=0.5 (Baseline)\nMean=%.2f, Std=%.2f', mean(max_streamflows_beta3), std(max_streamflows_beta3)));
xlabel('Maximum Streamflow (mm)');
ylabel('Probability');
grid on;

%% Part (d): Repeat with β = 7
fprintf('\nPart (d): Running simulations with β = 7...\n');

% Generate new rainfall sequences with β = 7
rainfall_data_beta7 = zeros(n_simulations, n_days);
beta_new = 7;

for sim = 1:n_simulations
    for day = 1:n_days
        if rand() < prob_rain
            rainfall_data_beta7(sim, day) = gamrnd(alpha, beta_new);
        else
            rainfall_data_beta7(sim, day) = 0;
        end
    end
end

max_streamflows_beta7 = zeros(1, n_simulations);

for sim = 1:n_simulations
    [streamflow, ~, ~, ~, ~] = simulate_streamflow(rainfall_data_beta7(sim,:), d, c, a, S_max, S_initial);
    max_streamflows_beta7(sim) = max(streamflow);
    
    if mod(sim, 200) == 0
        fprintf('Completed %d simulations for β=7\n', sim);
    end
end

fprintf('Mean maximum streamflow (β=7): %.2f mm\n', mean(max_streamflows_beta7));
fprintf('Std maximum streamflow (β=7): %.2f mm\n', std(max_streamflows_beta7));

% Plot for β=7
subplot(1,3,2);
histogram(max_streamflows_beta7, 30, 'Normalization', 'probability', 'FaceColor', 'red', 'FaceAlpha', 0.7);
title(sprintf('β=7, c=0.5\nMean=%.2f, Std=%.2f', mean(max_streamflows_beta7), std(max_streamflows_beta7)));
xlabel('Maximum Streamflow (mm)');
ylabel('Probability');
grid on;

%% Part (e): Repeat with c = 0.1
fprintf('\nPart (e): Running simulations with c = 0.1...\n');

c_new = 0.1;
max_streamflows_c01 = zeros(1, n_simulations);

for sim = 1:n_simulations
    [streamflow, ~, ~, ~, ~] = simulate_streamflow(rainfall_data(sim,:), d, c_new, a, S_max, S_initial);
    max_streamflows_c01(sim) = max(streamflow);
    
    if mod(sim, 200) == 0
        fprintf('Completed %d simulations for c=0.1\n', sim);
    end
end

fprintf('Mean maximum streamflow (c=0.1): %.2f mm\n', mean(max_streamflows_c01));
fprintf('Std maximum streamflow (c=0.1): %.2f mm\n', std(max_streamflows_c01));

% Plot for c=0.1
subplot(1,3,3);
histogram(max_streamflows_c01, 30, 'Normalization', 'probability', 'FaceColor', 'green', 'FaceAlpha', 0.7);
title(sprintf('β=3, c=0.1\nMean=%.2f, Std=%.2f', mean(max_streamflows_c01), std(max_streamflows_c01)));
xlabel('Maximum Streamflow (mm)');
ylabel('Probability');
grid on;

sgtitle('Comparison of Maximum Daily Streamflow Distributions');

%% Summary and Analysis
fprintf('\n=== SUMMARY STATISTICS ===\n');
fprintf('Case                    | Mean (mm) | Std (mm) | Min (mm) | Max (mm)\n');
fprintf('β=3, c=0.5 (baseline)  |  %6.2f  |  %5.2f  |  %5.2f  |  %5.2f\n', ...
    mean(max_streamflows_beta3), std(max_streamflows_beta3), min(max_streamflows_beta3), max(max_streamflows_beta3));
fprintf('β=7, c=0.5             |  %6.2f  |  %5.2f  |  %5.2f  |  %5.2f\n', ...
    mean(max_streamflows_beta7), std(max_streamflows_beta7), min(max_streamflows_beta7), max(max_streamflows_beta7));
fprintf('β=3, c=0.1             |  %6.2f  |  %5.2f  |  %5.2f  |  %5.2f\n', ...
    mean(max_streamflows_c01), std(max_streamflows_c01), min(max_streamflows_c01), max(max_streamflows_c01));

%% Analysis and Discussion
fprintf('\n=== ANALYSIS ===\n');
fprintf('Effect of changing β from 3 to 7:\n');
fprintf('- Increases mean rainfall intensity from %.1f to %.1f mm\n', alpha*beta, alpha*beta_new);
fprintf('- Increases variability in maximum streamflow\n');
fprintf('- Higher β leads to more extreme rainfall events\n\n');

fprintf('Effect of changing c from 0.5 to 0.1:\n');
fprintf('- Dramatically reduces maximum streamflow (%.1f%% decrease)\n', ...
    100*(1 - mean(max_streamflows_c01)/mean(max_streamflows_beta3)));
fprintf('- Lower baseflow coefficient means less drainage from storage\n');
fprintf('- Results in much lower overall streamflow generation\n');

% Save results to file
save('streamflow_results.mat', 'max_streamflows_beta3', 'max_streamflows_beta7', 'max_streamflows_c01', ...
     'rainfall_data', 'rainfall_data_beta7');

fprintf('\nResults saved to streamflow_results.mat\n');
fprintf('Analysis complete!\n');