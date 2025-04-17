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

to_physical = @(x_norm) lb + x_norm .* (ub - lb);
to_normalized = @(x_phys) (x_phys - lb) ./ (ub - lb);

%% Parameters algorithm
alpha = 5e-3;           % learning rate
penalty = 1e20;          % penalty multiplier
tol = 1e-3;             % tolerance for convergence
max_iter = 1000;         % max number of iterations

x0_norm = to_normalized(x0);  % safer initial guess
% Store the points for visualization
iter_points = []; 
%iter_points(end+1, :) = x0;
iter_points(end+1, :) = to_physical(x0_norm);  % store physical point

for iter = 1:max_iter
    x_phys = to_physical(x0_norm);
    
    f = objective(x_phys, W_base, rho);
    [c, ~] = nonlcon(x_phys, W_base, E, L, sigma_allow, disp_limit, ...
                     F_ref, node_coords, members, safety_fac);
    penalty_term = sum((max(0,c)).^2);
    F = f + penalty * penalty_term;

    grad = zeros(1, length(x0_norm));
    h = 1e-2;
    for i = 1:length(x0_norm)
        x_temp = x0_norm;
        x_temp(i) = x_temp(i) + h;
        x_phys_temp = to_physical(x_temp);

        f_temp = objective(x_phys_temp, W_base, rho);
        [c_temp, ~] = nonlcon(x_phys_temp, W_base, E, L, sigma_allow, disp_limit, ...
                              F_ref, node_coords, members, safety_fac);
        penalty_temp = sum((max(0,c_temp)).^2);
        F_temp = f_temp + penalty * penalty_temp;

        grad(i) = (F_temp - F) / h;
    end

    grad = grad / max(norm(grad), 1);  % normalize to avoid zigzagging
    x_new_norm = x0_norm - alpha * grad;

    % clamp to [0,1]
    x_new_norm = max(min(x_new_norm, 1), 0);

    % check convergence in physical space
    if norm(to_physical(x_new_norm) - to_physical(x0_norm)) < tol
        break;
    end

    iter_points(end+1, :) = to_physical(x_new_norm);  % store actual physical values
    x0_norm = x_new_norm;
end


x = to_physical(x0_norm);
mass_final = objective(x, W_base, rho);
fprintf('Optimal thickness t: %.4f m\n', x(1));
fprintf('Optimal height/width ratio r: %.4f\n', x(2));
fprintf('Minimum mass: %.4f kg\n', mass_final);


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

