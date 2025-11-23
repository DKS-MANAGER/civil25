% Basic MATLAB functionality test
try
    % Step 1: Clear workspace and setup
    clear; clc;
    disp('Step 1: Basic setup - OK');
    
    % Step 2: Test basic MATLAB functions
    x = 1:10;
    y = x.^2;
    disp('Step 2: Basic operations - OK');
    
    % Step 3: Test random number generation
    r = rand(1,10);
    disp('Step 3: Random numbers - OK');
    
    % Step 4: Test gamma distribution (requires Statistics Toolbox)
    g = gamrnd(1,1,1,10);
    disp('Step 4: Gamma distribution - OK');
    
    % Step 5: Test plotting
    figure('Visible', 'off');  % Create figure without displaying
    plot(x,y);
    close;
    disp('Step 5: Plotting - OK');
    
    disp('All basic tests passed successfully!');
    
catch ME
    disp('Error occurred in test:');
    disp(['Error message: ' ME.message]);
    disp(['Error occurred in: ' ME.stack(1).name]);
    disp(['Line number: ' num2str(ME.stack(1).line)]);
end