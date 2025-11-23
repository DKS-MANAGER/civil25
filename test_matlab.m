% BEGINNER'S MATLAB CODE FOR STREAMFLOW MODELING
% This code is written in a simple way for beginners to understand
% Each step is explained clearly with comments

% Clean up the workspace (like clearing a table before starting work)
clear all          % Remove all variables from memory
close all          % Close all figure windows
clc               % Clear the command window (makes it neat)

% STEP 1: SET UP THE PROBLEM
% These are our given values - like ingredients for cooking
number_of_simulations = 1000;    % How many times we run the experiment
number_of_days = 120;            % How many days in each simulation
chance_of_rain = 0.4;            % 40% chance of rain each day

% Gamma distribution parameters (describes how much it rains)
alpha = 5;                       % Shape parameter
beta = 3;                        % Scale parameter

% Model parameters (how the catchment behaves)
d = 0.3;                         % Surface runoff coefficient
c = 0.5;                         % Baseflow coefficient  
a = 0.6;                         % Infiltration coefficient
max_storage = 500;               % Maximum water storage (mm)
initial_storage = 100;           % Starting water storage (mm)

% Set random number seed (so results are the same each time)
rng(42);

% STEP 2: GENERATE RAINFALL DATA
% Create empty matrix to store rainfall (like empty boxes)
rainfall = zeros(number_of_simulations, number_of_days);

fprintf('Step 1: Generating rainfall data...\n');

% Loop through each simulation
for sim = 1:number_of_simulations
    
    % Loop through each day
    for day = 1:number_of_days
        
        % Flip a coin to see if it rains (random number between 0 and 1)
        random_number = rand();
        
        if random_number < chance_of_rain
            % It rains! Generate rainfall amount using gamma distribution
            rainfall_amount = gamrnd(alpha, beta);
            rainfall(sim, day) = rainfall_amount;
        else
            % No rain today
            rainfall(sim, day) = 0;
            end
        
        end
        end
        
    end
    
    % Show progress every 100 simulations
    if mod(sim, 100) == 0
        fprintf('Completed %d simulations\n', sim);
    end
    
end

% Check if our rainfall looks right
rainy_days = rainfall > 0;                    % Find all rainy days
actual_rain_chance = mean(rainy_days(:));     % Calculate actual probability
rainy_amounts = rainfall(rainfall > 0);       % Get only rainy day amounts

fprintf('\nChecking rainfall statistics:\n');
fprintf('Expected rain probability: 0.4, Actual: %.3f\n', actual_rain_chance);
fprintf('Expected rain amount mean: %.1f mm, Actual: %.1f mm\n', alpha*beta, mean(rainy_amounts));

% STEP 3: STREAMFLOW CALCULATION
fprintf('\nStep 2: Calculating streamflow for each simulation...\n');

% Create empty array to store maximum streamflow from each simulation
max_streamflow_beta3 = zeros(1, number_of_simulations);

% Loop through each simulation
for sim = 1:number_of_simulations
    
    % Get rainfall for this simulation
    daily_rainfall = rainfall(sim, :);
    
    % Create empty arrays for this simulation
    storage = zeros(1, number_of_days + 1);      % Need extra space for day 0
    surface_runoff = zeros(1, number_of_days);
    baseflow = zeros(1, number_of_days);
    infiltration = zeros(1, number_of_days);
    streamflow = zeros(1, number_of_days);
    
    % Set starting storage
    storage(1) = initial_storage;
    
    % Calculate streamflow for each day
    for day = 1:number_of_days
        
        % Calculate surface runoff (water that flows immediately)
        surface_runoff(day) = d * (storage(day) / max_storage) * daily_rainfall(day);
        
        % Calculate baseflow (water slowly draining from storage)
        baseflow(day) = c * storage(day);
        
        % Calculate infiltration (water soaking into ground)
        infiltration(day) = a * daily_rainfall(day);
        
        % Total streamflow = surface runoff + baseflow
        streamflow(day) = surface_runoff(day) + baseflow(day);
        
        % Update water storage for next day
        new_storage = storage(day) + infiltration(day) - baseflow(day);
        
        % Make sure storage doesn't go below 0 or above maximum
        if new_storage < 0
            storage(day + 1) = 0;
        elseif new_storage > max_storage
            storage(day + 1) = max_storage;
        else
            storage(day + 1) = new_storage;
        end
        
    end
    
    % Find the highest streamflow in this simulation
    max_streamflow_beta3(sim) = max(streamflow);
    
    % Show progress every 200 simulations
    if mod(sim, 200) == 0
        fprintf('Completed streamflow calculations for %d simulations\n', sim);
    end
    
end

% STEP 4: ANALYZE RESULTS FOR PART C
fprintf('\n=== PART C RESULTS ===\n');
mean_max_flow_beta3 = mean(max_streamflow_beta3);
std_max_flow_beta3 = std(max_streamflow_beta3);
min_max_flow_beta3 = min(max_streamflow_beta3);
max_max_flow_beta3 = max(max_streamflow_beta3);

fprintf('Statistics for maximum streamflow (beta=3, c=0.5):\n');
fprintf('Mean: %.2f mm\n', mean_max_flow_beta3);
fprintf('Standard deviation: %.2f mm\n', std_max_flow_beta3);
fprintf('Minimum: %.2f mm\n', min_max_flow_beta3);
fprintf('Maximum: %.2f mm\n', max_max_flow_beta3);

% Create histogram for part C
figure(1);
histogram(max_streamflow_beta3, 30, 'FaceColor', 'blue');
title('Distribution of Maximum Streamflow (beta=3, c=0.5)');
xlabel('Maximum Streamflow (mm)');
ylabel('Number of Simulations');
grid on;

% (Add similar steps for part d and e if needed)
