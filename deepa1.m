% ---- Read CSV file ----
filename = "deepa_grib_daily.csv";
data = readtable(filename);

% Convert time column from dd-MM-yyyy format
data.time = datetime(data.time, 'InputFormat', 'dd-MM-yyyy');

% ---- Select one date ----
selData = data(data.time == datetime(2024,1,15), :);

% Extract pressure & temperature
P_levels = selData.isobaricInhPa;       % hPa
T_levels = selData.temperature;         % Kelvin

% ---- Sort pressure in decreasing order ----
[~, idx] = sort(P_levels, 'descend');
P_levels = P_levels(idx);
T_levels = T_levels(idx);

% ---- Theory constants ----
g  = 9.8;       % m/s^2
Ra = 287;       % J/kg/K
alpha = 1;      % assumption (constant for now)
theoretical_exp = g / (alpha * Ra);

% ---- Compute β for each consecutive pair ----
nLayers = length(P_levels) - 1;

Results = table('Size',[nLayers 8], ...
                'VariableTypes',{'double','double','double','double','double','double','double','double'}, ...
                'VariableNames',{'P1','P2','T1','T2','Beta','Alpha','g_over_aRa','Difference'});

for i = 1:nLayers
    P1 = P_levels(i);
    P2 = P_levels(i+1);
    T1 = T_levels(i);
    T2 = T_levels(i+1);

    % Compute beta
    Beta = log(P2/P1) / log(T2/T1);

    % Difference
    diffVal = Beta - theoretical_exp;

    % Store in table (now includes Alpha)
    Results(i,:) = {P1,P2,T1,T2,Beta,alpha,theoretical_exp,diffVal};
end

% ---- Display results ----
disp(Results);

% ---- Plot with dual y-axes ----
figure;

% Left Y-axis for Beta and theoretical exponent
yyaxis left
plot(1:nLayers, Results.Beta, '-o', 'Color',[0 0 1], ...
     'MarkerFaceColor','b','MarkerSize',8,'LineWidth',2.5);  % No legend
hold on;
plot([1 nLayers], [theoretical_exp theoretical_exp], '-r', ...
     'LineWidth',2, 'DisplayName','g/(\alphaR_a)');          % Only ONE line drawn
ylabel('\beta and g/(\alphaR_a)');

% Right Y-axis for Difference
yyaxis right
plot(1:nLayers, Results.Difference, '-s', 'Color',[0 0.6 0], ...
     'MarkerFaceColor','g','MarkerSize',7,'LineWidth',2);    % No legend
ylabel('Difference (\beta - g/(\alphaR_a))');

% Labels and legend
xlabel('Layer index (1 = top layer)');
title('Comparison of \beta, \alpha, Theoretical Exponent, and Difference');
grid on;
legend('show','Location','southoutside','Orientation','horizontal');
hold off;