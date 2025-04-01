clc; clear; close all;


%% using a Gradient-Based Constrained Nonlinear Optimization Algorithm (wint penalization)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Nodes and Members
node_coords = [...
     0    0;
     7    0;
    14    0;
    21    0;
    3.5   6;
   10.5   6;
   17.5   6];

members = [...
    1 2;
    2 3;
    3 4;
    5 6;
    6 7;
    1 5;
    2 5;
    2 6;
    3 6;
    3 7;
    4 7];


%% Properties
rho = 7850;        % kg/m^3
E = 210e9;         % Pa
sigma_allow = 250e6; % Pa
W_base = 0.3;      % m (constant width)
disp_limit = 0.02; % m (20 mm max displacement)




%% Plot b ridge
plot_bridge;






%% Optimization Setup
x0 = [0.05, 2.5];      % Initial guess [t, r]
lb = [0.01, 1];       % Lower bounds [t_min, r_min]
ub = [0.5, 10];        % Upper bounds [t_max, r_max]


%% Parameters algorithm
alpha = 1e-3;           % learning rate
penalty = 1e10;          % penalty multiplier
tol = 1e-9;             % tolerance for convergence
max_iter = 500;         % max number of iterations

for iter = 1:max_iter
    % obj and constraints
  
    f = objective(x0,W_base,rho);
    [c,~] = nonlcon(x0, W_base, E, sigma_allow, disp_limit,node_coords,members);

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
        [c_temp, ~] = nonlcon(x_temp, W_base, E, sigma_allow, disp_limit, node_coords, members);
        F_temp = f_temp + penalty * sum((max(0, c_temp)).^2);

        grad(i) = (F_temp - F) / h;
    end
    
    % gradient descent step
    x_new = x0 - alpha*grad;

    %project into bounds
    x_new = max(min(x_new,ub), lb);

    %check convergence
    if norm(x_new - x0) < tol
        break;
    end
    
    x = x_new;

end

%% CHeck if constraints are satisfied

[c_final, ~] = nonlcon(x, W_base, E, sigma_allow, disp_limit, node_coords, members);
disp('Constraint values at optimum (should all be ≤ 0):');
disp(c_final);




%% Final result
mass_final = objective(x, W_base, rho);
disp('--- GRADIENT DESCENT RESULT ---');
fprintf('Optimal thickness t: %.4f m\n', x(1));
fprintf('Optimal height/width ratio r: %.4f\n', x(2));
fprintf('Minimum mass: %.4f kg\n', mass_final);



%% Final Contour Plot
% figure;
% hold on;
% 
% % Plot mass contours
% contourf(T, R, MASS, 30, 'LineColor', 'none');
% colormap('parula');
% colorbar;
% xlabel('Thickness t (m)');
% ylabel('Height-to-width ratio r');
% title('Mass Contour with Optimization Path');
% 
% % Constraint boundaries overlay
% contour(T, R, CONSTRAINT_VIOLATION, [1 1], 'k--', 'LineWidth', 2);
% 
% % If you store path, plot it here
% if exist('opt_path', 'var') && ~isempty(opt_path)
%     plot(opt_path(:,1), opt_path(:,2), 'w.-', 'LineWidth', 2, 'MarkerSize', 10);
% end
% 
% % Plot final point
% plot(x(1), x(2), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
% text(x(1) + 0.01, x(2), 'Optimal Point', 'FontSize', 10, 'Color', 'r');
% 
% hold off;
