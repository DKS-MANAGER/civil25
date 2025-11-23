% Filename: compute_lapse_rate.m
% Description: Compute daily bulk lapse rate from ERA5 data in Excel.

% Specify Excel file
filename = 'era5_corrected_complete.xlsx';  % Replace with actual filename

% Read the data
% Assumes sheet has columns: Day, Pressure, Temperature, GeoHeight
T = readtable(filename);

% Unique days in the dataset
days = unique(T.Day);

% Preallocate results
n = numel(days);
bulkLapseRate = nan(n,1);

% Loop over each day
for i = 1:n
    day = days(i);
    
    % Extract data for this day
    idx = T.Day == day;
    data = T(idx, :);
    
    % Identify surface (max Pressure) and top (min Pressure)
    [~, surfIdx] = max(data.Pressure);
    [~, topIdx]  = min(data.Pressure);
    
    T_surf = data.Temperature(surfIdx);
    T_top  = data.Temperature(topIdx);
    Z_surf = data.GeoHeight(surfIdx);
    Z_top  = data.GeoHeight(topIdx);
    
    % Compute bulk lapse rate (K/m)
    bulkLapseRate(i) = - (T_top - T_surf) / (Z_top - Z_surf);
end

% Combine results into a table
results = table(days, bulkLapseRate, 'VariableNames', {'Day','BulkLapseRate_K_per_m'});

% Display results
disp(results);

% Optionally, write results back to Excel
writetable(results, 'daily_lapse_rates.xlsx', 'Sheet',1);
