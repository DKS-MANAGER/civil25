function [results, beta_theoretical] = beta_analysis()
% BETA_ANALYSIS Complete Beta Exponent Analysis Function
% This function performs comprehensive beta exponent analysis for 
% atmospheric pressure-temperature relationships using Ra = rd(1 + 0.61q)
%
% Outputs:
%   results - Table with all pressure pair analysis results
%   beta_theoretical - Theoretical beta value using moist air constant
%
% Example usage:
%   [results, beta_theo] = beta_analysis();
%   fprintf('Theoretical beta: %.3f\n', beta_theo);

% Constants
g = 9.81;           % acceleration due to gravity (m/s²)
alpha = 0.0065;     % standard atmospheric lapse rate (K/m)
rd = 287.04;        % specific gas constant of dry air (J/(kg·K))

fprintf('🔬 Complete Beta Exponent Analysis Solution\n');
fprintf('Using corrected Ra values with moisture effects\n\n');

% Load data (assuming CSV format for MATLAB)
try
    temp_data = readtable('era5_corrected_complete.csv');
    q_data = readtable('qvalue.csv');
    fprintf('✅ Data loaded successfully\n');
catch ME
    fprintf('❌ Error loading data: %s\n', ME.message);
    return;
end

% Calculate Ra for moist air
q_data.Ra = rd * (1 + 0.61 * q_data.SpecificHumidity);

% Convert date strings to datetime if needed
temp_data.Date = datetime(temp_data.Day);
q_data.Date = datetime(q_data.Day);

% Merge datasets
merged_data = innerjoin(temp_data, q_data(:,{'Date','Pressure','Ra'}), 'Keys', {'Date','Pressure'});

% Calculate theoretical beta
mean_Ra = mean(merged_data.Ra);
beta_theoretical = g / (alpha * mean_Ra);

fprintf('Average Ra (moist air): %.1f J/(kg·K)\n', mean_Ra);
fprintf('Theoretical β = %.3f\n\n', beta_theoretical);

% Filter data for 15th of each month
day_15_mask = day(merged_data.Date) == 15;
day_15_data = merged_data(day_15_mask, :);

% Get unique pressure levels
pressure_levels = unique(day_15_data.Pressure);
pressure_levels = sort(pressure_levels, 'descend'); % Sort high to low

% Initialize results arrays
P1_list = [];
P2_list = [];
P_ratio_list = [];
T_ratio_list = [];
beta_calc_list = [];
beta_theo_list = [];
diff_list = [];
region_list = {};

% Calculate beta exponents for all pressure pairs
fprintf('Calculating beta exponents for pressure pairs...\n');
pair_count = 0;

for i = 1:length(pressure_levels)
    for j = 1:length(pressure_levels)
        p1 = pressure_levels(i);
        p2 = pressure_levels(j);

        if p1 > p2  % Only consider pairs where P1 > P2
            % Get data for each pressure level
            data_p1 = day_15_data(day_15_data.Pressure == p1, :);
            data_p2 = day_15_data(day_15_data.Pressure == p2, :);

            % Find common dates
            [common_dates, idx1, idx2] = intersect(data_p1.Date, data_p2.Date);

            if ~isempty(common_dates)
                % Calculate ratios
                P_ratio = p1 / p2;
                T_ratio = mean(data_p1.Temperature(idx1)) / mean(data_p2.Temperature(idx2));
                Ra_avg = (mean(data_p1.Ra(idx1)) + mean(data_p2.Ra(idx2))) / 2;

                if T_ratio > 1
                    % Calculate beta exponent
                    beta_calculated = log(P_ratio) / log(T_ratio);
                    beta_theoretical_pair = g / (alpha * Ra_avg);
                    difference = beta_calculated - beta_theoretical_pair;

                    % Determine atmospheric region
                    if p2 >= 200
                        region = 'Troposphere';
                    elseif p2 >= 100
                        region = 'Tropopause';
                    else
                        region = 'Stratosphere';
                    end

                    % Store results
                    pair_count = pair_count + 1;
                    P1_list(pair_count) = p1;
                    P2_list(pair_count) = p2;
                    P_ratio_list(pair_count) = P_ratio;
                    T_ratio_list(pair_count) = T_ratio;
                    beta_calc_list(pair_count) = beta_calculated;
                    beta_theo_list(pair_count) = beta_theoretical_pair;
                    diff_list(pair_count) = difference;
                    region_list{pair_count} = region;
                end
            end
        end
    end
end

% Create results table
results = table(P1_list', P2_list', P_ratio_list', T_ratio_list', ...
                beta_calc_list', beta_theo_list', diff_list', region_list', ...
                'VariableNames', {'P1_hPa', 'P2_hPa', 'P1_P2_ratio', 'T1_T2_ratio', ...
                                 'Beta_calculated', 'Beta_theoretical', 'Difference', 'Region'});

% Find best agreements
[~, sort_idx] = sort(abs(results.Difference));
best_results = results(sort_idx(1:10), :);

fprintf('\n🏆 Top 10 Best Beta Exponent Agreements:\n');
fprintf('Pair\t\tβ_calc\tβ_theo\tDiff\n');
fprintf('----\t\t------\t------\t----\n');
for i = 1:height(best_results)
    fprintf('%d-%d\t\t%.3f\t%.3f\t%.3f\n', ...
            best_results.P1_hPa(i), best_results.P2_hPa(i), ...
            best_results.Beta_calculated(i), best_results.Beta_theoretical(i), ...
            best_results.Difference(i));
end

% Save results to CSV
writetable(results, 'beta_analysis_matlab_results.csv');
fprintf('\n✅ Results saved to beta_analysis_matlab_results.csv\n');

% Summary statistics by region
fprintf('\n📋 Analysis by Atmospheric Region:\n');
regions = unique(results.Region);
for i = 1:length(regions)
    region_mask = strcmp(results.Region, regions{i});
    region_data = results(region_mask, :);

    mean_beta = mean(region_data.Beta_calculated);
    mean_diff = mean(region_data.Difference);
    count = sum(region_mask);

    status = '✅';
    if abs(mean_diff) >= 1.0
        status = '❌';
    end

    fprintf('%s: %d pairs, Mean β = %.3f, Mean diff = %.3f %s\n', ...
            regions{i}, count, mean_beta, mean_diff, status);
end

fprintf('\n🎯 Analysis complete!\n');

end

% Helper function for data visualization (optional)
function plot_beta_results(results)
    % Create scatter plot of beta values
    figure('Name', 'Beta Exponent Analysis', 'Position', [100, 100, 1200, 800]);

    subplot(2, 2, 1);
    scatter(1:height(results), results.Beta_calculated, 50, 'filled');
    hold on;
    yline(mean(results.Beta_theoretical), 'r--', 'LineWidth', 2);
    xlabel('Pressure Pair Index');
    ylabel('Beta Exponent');
    title('Calculated Beta Exponents');
    grid on;

    subplot(2, 2, 2);
    histogram(results.Difference, 20);
    xlabel('Difference (Calc - Theo)');
    ylabel('Frequency');
    title('Distribution of Differences');
    grid on;

    subplot(2, 2, 3);
    scatter(results.P1_P2_ratio, results.T1_T2_ratio, 50, results.Beta_calculated, 'filled');
    xlabel('Pressure Ratio (P1/P2)');
    ylabel('Temperature Ratio (T1/T2)');
    title('Pressure vs Temperature Ratios');
    colorbar;
    grid on;

    subplot(2, 2, 4);
    boxplot(results.Beta_calculated, results.Region);
    ylabel('Beta Exponent');
    title('Beta by Atmospheric Region');
    grid on;

    sgtitle('Complete Beta Exponent Analysis Results');
end