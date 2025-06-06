clc; clear; close all;


%% using a Gradient-Based Constrained Nonlinear Optimization Algorithm (wint penalization)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Properties and Nodes
params;


%% Plot b ridge
plot_bridge;




%x_slp = [0.4, 1.0];   % Initial guess [t, r]
x_slp = [0.06, 8];


lb = [0.005, 1];       % Lower bounds [t_min, r_min]
ub = [0.03, 10];        % Upper bounds [t_max, r_max]
% thickness between 5 and 30 mm
% r = H/W

trust_region = 0.1;   % Max allowed change per iteration
max_iter = 50;
tol =  1e-3;

% Store the points for visualization
iter_points = []; 
iter_points(end+1, :) = x_slp;

for iter = 1:max_iter
    %  objective and constraints linearization

    f0 = objective(x_slp, W_base, rho); %objective at current x
    grad_f = finite_diff(@(x) objective(x, W_base, rho), x_slp); 

    [c, ~] = nonlcon(x_slp, W_base, E, L, sigma_allow, disp_limit, F_ref,node_coords, members,safety_fac);
    A = zeros(length(c), length(x_slp)); % jacobian constraints
    for j = 1:length(c)
        A(j,:) = finite_diff(@(x) nonlcon(x, W_base, E, L, sigma_allow, disp_limit, F_ref,node_coords, members,safety_fac), x_slp, j);
    end

    % LP setup: minimize grad_f * dx
    % subject to: A*dx + c <= 0, bounds, and trust region

     % trust region formulation
    dx_lb = max(lb - x_slp, -trust_region);
    dx_ub = min(ub - x_slp, trust_region);

    % Enforce that lower bound is not greater than upper bound
    invalid_bounds = dx_lb > dx_ub;
    dx_lb(invalid_bounds) = dx_ub(invalid_bounds);  % set dx_lb = dx_ub if invalid    

    fprintf('x = [%f, %f]\n', x_slp(1), x_slp(2));
    fprintf('dx_lb = [%f, %f], dx_ub = [%f, %f]\n', dx_lb(1), dx_lb(2), dx_ub(1), dx_ub(2));

    % Solve LP using linprog
    fprintf('lb = [%.4f, %.4f], ub = [%.4f, %.4f]\n', lb(1), lb(2), ub(1), ub(2));
    options = optimoptions('linprog','Display','iter');
    [dx, fval, exitflag] = linprog(grad_f, A, -c, [], [], dx_lb, dx_ub, options);

    if exitflag ~= 1
        fprintf('LP failed at iter %d\n', iter);
        break;
    end

    % feasibility check
    dx = max(min(dx, dx_ub), dx_lb);

    x_new = x_slp + dx(1,:);
    x_new = max(min(x_new, ub), lb); 
    [c_new, ~] = nonlcon(x_new, W_base, E, L, sigma_allow, disp_limit, F_ref,node_coords, members,safety_fac);

    if any(c_new > 0)
        fprintf('⚠️ Step violates constraints. Reducing trust region.\n');
        trust_region = trust_region * 0.5;
        continue; % Retry iteration with smaller trust region
    end



    fprintf('Iter %d → x = [%.4f, %.4f], f = %.2f\n', iter, x_new(1), x_new(2), objective(x_new, W_base, rho));
    % convergence check

    if norm(x_new - x_slp) < tol
        % fprintf('Step size: %.10f\n', norm(x_new - x_slp));
        break;
    end
    iter_points(end+1, :) = x_new;
    x_slp = x_new;
end

% final results
mass_slp = objective(x_slp, W_base, rho);
[c_slp, ~] = nonlcon(x_slp, W_base, E, L, sigma_allow, disp_limit, F_ref,node_coords, members,safety_fac);

disp('--- SLP RESULT ---');
fprintf('Optimal thickness t: %.4f m\n', x_slp(1));
fprintf('Optimal height/width ratio r: %.4f\n', x_slp(2));
fprintf('Minimum mass: %.4f kg\n', mass_slp);
disp('Constraint values at optimum (should all be ≤ 0):');
disp(c_slp);

%Plot points
figure(2)
plot_optimization_contour(iter_points)

figure(3)
visualize_linearized_constraints(x_slp)



