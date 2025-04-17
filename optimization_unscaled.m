clc; clear; close all;


%% using a Gradient-Based Constrained Nonlinear Optimization Algorithm (wint penalization)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Properties
params;


%% Plot b ridge
plot_bridge;


%% Optimization Setup
x0 = [0.06, 8];      % Initial guess [t, r]
lb = [0.005, 1];       % Lower bounds [t_min, r_min]
ub = [0.03, 10];        % Upper bounds [t_max, r_max]
% thickness between 5 and 30 mm
% r = H/W


%% Parameters algorithm
alpha = 1e-3;           % learning rate
penalty = 1e20;          % penalty multiplier
tol =  1e-3;             % tolerance for convergence
max_iter = 1000;         % max number of iterations


% Store the points for visualization
iter_points = []; 
iter_points(end+1, :) = x0;



for iter = 1:max_iter
    % obj and constraints
    alpha = 1e-3 / sqrt(iter);
    f = objective(x0,W_base,rho);
    [c,~] = nonlcon(x0, W_base, E,L, sigma_allow, disp_limit,F_ref,node_coords,members,safety_fac);

    % penalty
    penalty_term  = sum((max(0,c)).^2);
    F = f + penalty* penalty_term;

    %difference gradient
    grad = zeros(1, length(x0));
    h = 1e-2;
    for i = 1:length(x0)
        x_temp = x0;
        x_temp(i) = x_temp(i) + h;

        f_temp = objective(x_temp, W_base, rho);
        [c_temp, ~] = nonlcon(x_temp, W_base, E,L, sigma_allow, disp_limit, F_ref,node_coords, members,safety_fac);
        F_temp = f_temp + penalty * sum((max(0, c_temp)).^2);

        grad(i) = (F_temp - F) / h;
    end

    %grad(1) = grad(1) / (ub(1) - lb(1));   % scale t gradient
    %grad(2) = grad(2) / (ub(2) - lb(2));   % scale r gradient


    % gradient descent step
    x_new = x0 - alpha*grad;

    %project into bounds
    x_new = max(min(x_new,ub), lb);

    %check convergence
    if norm(x_new - x0) < tol
        break;
    end
    iter_points(end+1, :) = x_new;
    x = x_new;
    x0 = x_new;

end


%% CHeck if constraints are satisfied

[c_final, ~] = nonlcon(x, W_base, E, L, sigma_allow, disp_limit, F_ref,node_coords, members,safety_fac);
disp('Constraint values at optimum (should all be ≤ 0):');
disp(c_final);




%% Final result
mass_final = objective(x, W_base, rho);
disp('--- GRADIENT DESCENT RESULT ---');
fprintf('Optimal thickness t: %.4f m\n', x(1));
fprintf('Optimal height/width ratio r: %.4f\n', x(2));
fprintf('Minimum mass: %.4f kg\n', mass_final);

%Plot points
plot_optimization_contour(iter_points)


