clc;
clearvars;

% Given data points
x = [0 1 2 3 4];
y = [1 2 1 3 2];

% Points at which interpolation is required
x_interp = [2.5 3.5];

%  Linear Interpolation 
fprintf('Linear Interpolation Results:\n');
for k = 1:length(x_interp)
    % Locate interval where x lies
    for i = 1:length(x)-1
        if x_interp(k) >= x(i) && x_interp(k) <= x(i+1)
            % Apply linear interpolation formula
            f_val = y(i) + ((x_interp(k) - x(i)) / (x(i+1) - x(i))) * (y(i+1) - y(i));
            fprintf('f(%.1f) = %.4f\n', x_interp(k), f_val);
            break;
        end
    end
end

% Lagrange Polynomial Interpolation 
fprintf('\nLagrange Interpolation Results:\n');
for k = 1:length(x_interp)
    f_lagr = 0;
    for i = 1:length(x)
        L = 1;
        for j = 1:length(x)
            if j ~= i
                L = L * (x_interp(k) - x(j)) / (x(i) - x(j));
            end
        end
        f_lagr = f_lagr + y(i) * L;
    end
    fprintf('f(%.1f) = %.4f\n', x_interp(k), f_lagr);
end

%  Visualization 
xx = linspace(min(x), max(x), 100);
yy_lagr = zeros(size(xx));

for k = 1:length(xx)
    f_temp = 0;
    for i = 1:length(x)
        L = 1;
        for j = 1:length(x)
            if j ~= i
                L = L * (xx(k) - x(j)) / (x(i) - x(j));
            end
        end
        f_temp = f_temp + y(i) * L;
    end
    yy_lagr(k) = f_temp;
end

% Plot results
plot(x, y, 'ro', 'MarkerFaceColor', 'r');
hold on;
plot(xx, yy_lagr, 'b-', 'LineWidth', 1.5);
plot(x_interp, interp1(x, y, x_interp), 'gs', 'MarkerFaceColor', 'g');
xlabel('x');
ylabel('f(x)');
title('Interpolation Comparison');
legend('Data points', 'Lagrange Polynomial', 'Linear Estimation', 'Location', 'best');
grid on;
