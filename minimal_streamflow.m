% Minimal Streamflow Model
function minimal_streamflow()
    try
        % 1. Setup
        clear; clc;
        rng(42);
        disp('1. Setup complete');
        
        % 2. Parameters
        params = setup_parameters();
        disp('2. Parameters initialized');
        
        % 3. Generate sample rainfall
        rainfall = generate_rainfall(params);
        disp('3. Rainfall generated');
        
        % 4. Calculate streamflow
        [flow, ~] = calculate_streamflow(rainfall, params);
        disp('4. Streamflow calculated');
        
        % 5. Plot results
        plot_results(rainfall, flow);
        disp('5. Results plotted');
        
        disp('Program completed successfully!');
        
    catch ME
        disp('Error occurred:');
        disp(['Message: ' ME.message]);
        disp(['Location: ' ME.stack(1).name ' line ' num2str(ME.stack(1).line)]);
    end
end

function params = setup_parameters()
    params.n_days = 10;  % Reduced for testing
    params.prob_rain = 0.4;
    params.alpha = 5;
    params.beta = 3;
    params.d = 0.3;
    params.c = 0.5;
    params.a = 0.6;
    params.S_max = 500;
    params.S_initial = 100;
end

function rainfall = generate_rainfall(params)
    rainfall = zeros(1, params.n_days);
    for day = 1:params.n_days
        if rand() < params.prob_rain
            % Use local gamma sampler to avoid dependency on Statistics Toolbox
            rainfall(day) = gamrnd_local(params.alpha, params.beta);
        end
    end
end

function [streamflow, storage] = calculate_streamflow(rainfall, params)
    n_days = length(rainfall);
    storage = zeros(1, n_days + 1);
    streamflow = zeros(1, n_days);
    storage(1) = params.S_initial;
    
    for t = 1:n_days
        % Surface runoff
        sr = params.d * (storage(t) / params.S_max) * rainfall(t);
        
        % Baseflow
        bf = params.c * storage(t);
        
        % Infiltration
        inf = params.a * rainfall(t);
        
        % Total streamflow
        streamflow(t) = sr + bf;
        
        % Update storage
        new_storage = storage(t) + inf - bf;
        storage(t + 1) = max(0, min(new_storage, params.S_max));
    end
    storage = storage(2:end);
end

function plot_results(rainfall, streamflow)
    figure('Position', [100, 100, 800, 400]);
    
    % Plot rainfall and streamflow
    days = 1:length(rainfall);
    
    subplot(2,1,1);
    bar(days, rainfall);
    title('Daily Rainfall');
    ylabel('Rainfall (mm)');
    grid on;
    
    subplot(2,1,2);
    plot(days, streamflow, 'b-', 'LineWidth', 1.5);
    title('Daily Streamflow');
    xlabel('Day');
    ylabel('Streamflow (mm)');
    grid on;
    
    % Save figure for headless runs
    try
        drawnow;
        saveas(gcf, 'minimal_streamflow_plot.png');
    catch
        % ignore save errors in restricted environments
    end
end

%% Local helper: Gamma random number generator (Marsaglia and Tsang)
function r = gamrnd_local(k, theta, varargin)
% GAMRND_LOCAL  Generate Gamma random numbers with shape k and scale theta
% Supports scalar or array output using optional size arguments, e.g. gamrnd_local(k,theta,m,n)
% This implementation uses the Marsaglia and Tsang method (2000) for k>0.

% Parse output size
if isempty(varargin)
    outSize = [1 1];
else
    outSize = cell2mat(varargin);
end

% Number of samples to generate
num = prod(outSize);

% Preallocate
r = zeros(outSize);

% Handle special case k==1 -> exponential
if k == 1
    r = theta * (-log(rand(outSize)));
    return;
end

% Use Marsaglia & Tsang method for k>0
for idx = 1:num
    kk = k;
    % For shapes < 1, use transformation
    if kk < 1
        % Gamma(k) = Gamma(k+1) * U^(1/k)
        % Generate Gamma(k+1) then multiply by U^(1/k)
        kk1 = kk + 1;
        % generate gamma with shape kk1
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