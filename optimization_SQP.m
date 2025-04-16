clc; clear; close all;

%% Load Parameters and Truss
params;
plot_bridge;

%% Initial Design and Bounds
x0 = [0.06, 8];            % Initial guess: [thickness t, height/width ratio r]
lb = [0.005, 1];              % Lower bounds
ub = [0.03, 10];              % Upper bounds

% Store the points for visualization
iter_points = []; 
iter_points(end+1, :) = x0;


%% Set SQP Optimizer Options
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'sqp', ...
    'MaxIterations', 100, ...
    'OptimalityTolerance', 1e-3, ...
    'StepTolerance', 1e-3, ...
    'ConstraintTolerance', 1e-3);

%% Run Optimization with fmincon
[x_opt, mass_opt, exitflag, output, lambda] = fmincon( ...
    @(x) objective(x, W_base, rho), ...                     % Objective function
    x0, [], [], [], [], lb, ub, ...                         % No linear constraints
    @(x) nonlcon(x, W_base, E, L, sigma_allow, disp_limit, ...
                 F_ref, node_coords, members, safety_fac), ... % Nonlinear constraints
    options);
iter_points(end+1, :) = x_opt;


%% Output Results
fprintf('\n--- SQP Optimization Results ---\n');
fprintf('Optimal thickness t: %.4f m\n', x_opt(1));
fprintf('Optimal height/width ratio r: %.4f\n', x_opt(2));
fprintf('Minimum mass: %.4f kg\n', mass_opt);

[c_opt, ~] = nonlcon(x_opt, W_base, E, L, sigma_allow, disp_limit, ...
                     F_ref, node_coords, members, safety_fac);

disp('Constraint values at optimum (should all be ≤ 0):');
disp(c_opt);

%% Optional: Plot Optimization Point on Contour
plot_optimization_contour(iter_points);
