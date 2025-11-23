% Simplified test of stochastic streamflow modeling
clear all; close all; clc;

try
    % Test 1: Basic setup and random number generation
    disp('Test 1: Setting up parameters...');
    rng(42);
    n_simulations = 10; % Reduced for testing
    n_days = 120;
    prob_rain = 0.4;
    alpha = 5;
    beta = 3;
    d = 0.3;
    c = 0.5;
    a = 0.6;
    S_max = 500;
    S_initial = 100;
    
    % Test 2: Generate small rainfall dataset
    disp('Test 2: Generating rainfall data...');
    rainfall_data = zeros(n_simulations, n_days);
    for sim = 1:n_simulations
        for day = 1:n_days
            if rand() < prob_rain
                rainfall_data(sim, day) = gamrnd_local(alpha, beta);
            end
        end
    end
    disp('Successfully generated rainfall data');
    
    % Test 3: Test streamflow calculation for one series
    disp('Test 3: Testing streamflow calculation...');
    [streamflow, storage, surface_runoff, baseflow, infiltration] = simulate_streamflow(rainfall_data(1,:), d, c, a, S_max, S_initial);
    disp(['Maximum streamflow: ' num2str(max(streamflow))]);
    
    % Test 4: Basic plotting
    disp('Test 4: Creating test plot...');
    figure;
    plot(1:n_days, streamflow);
    title('Test Streamflow Series');
    xlabel('Day');
    ylabel('Streamflow (mm)');
    
    disp('All tests completed successfully!');
    
catch ME
    disp(['Error occurred: ' ME.message]);
    disp(['in ' ME.stack(1).name ' at line ' num2str(ME.stack(1).line)]);
end

% Function definition at the end
function [streamflow, storage, surface_runoff, baseflow, infiltration] = simulate_streamflow(rainfall_series, d, c, a, S_max, S_initial)
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
        % Calculate surface runoff
        surface_runoff(t) = d * (storage(t) / S_max) * rainfall_series(t);
        
        % Calculate baseflow
        baseflow(t) = c * storage(t);
        
        % Calculate infiltration
        infiltration(t) = a * rainfall_series(t);
        
        % Calculate streamflow
        streamflow(t) = surface_runoff(t) + baseflow(t);
        
        % Update storage
        new_storage = storage(t) + infiltration(t) - baseflow(t);
        storage(t + 1) = max(0, min(new_storage, S_max));
    end
    
    % Return only the storage values for days 1 to n_days
    storage = storage(2:end);
end