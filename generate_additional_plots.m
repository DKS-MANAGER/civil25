% Generate additional plots (CDF and return-period) from streamflow_results.mat
try
    load('streamflow_results.mat', 'max_streamflows_beta3', 'max_streamflows_beta7', 'max_streamflows_c01');
catch
    error('streamflow_results.mat not found. Run sto4a_fixed.m first.');
end

% Ensure column vectors
b3 = max_streamflows_beta3(:);
b7 = max_streamflows_beta7(:);
c01 = max_streamflows_c01(:);

% CDF plot (manual empirical CDF to avoid Statistics Toolbox)
% CDF plot (separate figure)
figure('Position',[100,100,600,400]);
[x3, idx3] = sort(b3);
f3 = (1:length(x3))'./length(x3);
[x7, idx7] = sort(b7);
f7 = (1:length(x7))'./length(x7);
[x01, idx01] = sort(c01);
f01 = (1:length(x01))'./length(x01);
plot(x3,f3,'b-','LineWidth',1.5); hold on;
plot(x7,f7,'r-','LineWidth',1.5);
plot(x01,f01,'g-','LineWidth',1.5);
legend({'\beta=3','\beta=7','c=0.1'},'Location','best');
xlabel('Maximum daily streamflow (mm)'); ylabel('Empirical CDF');
title('Empirical CDF of maximum daily streamflow'); grid on;
try
    saveas(gcf,'sto4a_cdf.png');
catch
    warning('Could not save sto4a_cdf.png');
end
close(gcf);

% Return period plot (separate figure)
figure('Position',[100,100,600,400]);
s3 = sort(b3,'descend');
prob = (1:length(s3))'./(length(s3)+1);
ret3 = 1./prob;
semilogy(ret3,s3,'b-','LineWidth',1.2); hold on;
s7 = sort(b7,'descend'); prob7 = (1:length(s7))'./(length(s7)+1); semilogy(1./prob7,s7,'r-','LineWidth',1.2);
s01 = sort(c01,'descend'); prob01 = (1:length(s01))'./(length(s01)+1); semilogy(1./prob01,s01,'g-','LineWidth',1.2);
legend({'\beta=3','\beta=7','c=0.1'},'Location','best');
xlabel('Return period (1/(1-F))'); ylabel('Maximum daily streamflow (mm)');
title('Return-period plot (log-scale)'); grid on;
try
    saveas(gcf,'sto4a_return_period.png');
catch
    warning('Could not save sto4a_return_period.png');
end
close(gcf);

fprintf('Additional plots saved: sto4a_cdf.png, sto4a_return_period.png\n');
