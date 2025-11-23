% Read the CSV file
filename = "C:\Users\DKS\Desktop\p1\satish_grib_daily.csv";
data = readtable(filename);

% Convert time column from dd-MM-yyyy format
data.time = datetime(data.time, 'InputFormat', 'dd-MM-yyyy');

% Add separate columns for year, month, and day
data.Year  = year(data.time);
data.Month = month(data.time);
data.Day   = day(data.time);

% ---- Filter data for the 15th day of each month ----
day15 = data(data.Day == 15 & data.Year == 2024, :);  % only 2024

% Unique months available
monthsAvailable = unique(day15.Month);

% Create figure
figure;
hold on;

% Loop through each month and plot
for i = 1:length(monthsAvailable)
    thisMonth = monthsAvailable(i);
    
    % Extract data for the 15th of this month
    subset = day15(day15.Month == thisMonth, :);
    
    % Sort by pressure level (descending)
    subset = sortrows(subset, "isobaricInhPa", "descend");
    
    % Plot pressure vs. geopotential height
    plot(subset.isobaricInhPa, subset.geopotential_height, '-o', ...
         'DisplayName', datestr(datetime(2024,thisMonth,15), 'mmm'));
end

% Reverse x-axis so pressure decreases left to right
set(gca, 'XDir','reverse');

% Labels and title
xlabel('Pressure Level (hPa)');
ylabel('Geopotential Height (m)');
title('Pressure vs. Geopotential Height on 15th of Each Month (2024)');

% Show legend at bottom
legend('show', 'Location', 'southoutside', 'Orientation', 'horizontal');

% Force full values on y-axis
ax = gca;
ax.YAxis.Exponent = 0;

grid on;