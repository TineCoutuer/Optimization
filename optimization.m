clc; clear; close all;


%% using a Gradient-Based Constrained Nonlinear Optimization Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Properties
rho = 7850;        % kg/m^3
E = 210e9;         % Pa
sigma_allow = 250e6; % Pa
W_base = 0.3;      % m (constant width)
disp_limit = 0.02; % m (20 mm max displacement)

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


%% Plot b ridge
plot_bridge;

%% Define Design Space for Contour Plot
t_vals = linspace(0.02, 0.3, 50);   % Thickness range (m)
r_vals = linspace(1.0, 5.0, 50);    % Height-to-width ratio range
[T, R] = meshgrid(t_vals, r_vals);
MASS = zeros(size(T));              % Objective function (mass)
CONSTRAINT_VIOLATION = zeros(size(T)); % Constraint map

for i = 1:numel(T)
    t = T(i);
    r = R(i);
    
    % Compute Mass
    W = W_base;
    H = r * W_base;
    A = 2*W*t + (H-2*t)*t;
    L_A = 7;
    L_B = sqrt((L_A/2)^2 + H^2);
    MASS(i) = rho * A * (5*L_A + 6*L_B);
    
    % Check Constraint Violations
    [c, ~] = nonlcon([t, r], W_base, E, sigma_allow, disp_limit, node_coords, members);
    
    if any(c > 0)
        CONSTRAINT_VIOLATION(i) = 1; % If violated
    else
        CONSTRAINT_VIOLATION(i) = 0;
    end
end


%% Optimization Setup
x0 = [0.05, 2.5];      % Initial guess [t, r]
lb = [0.005, 1];       % Lower bounds [t_min, r_min]
ub = [0.5, 10];        % Upper bounds [t_max, r_max]


options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
    'OutputFcn', @optimize_plot_callback); % Call output function to store iterations

global opt_path;
opt_path = []; % Store optimization path

[x_opt, mass_opt] = fmincon(@(x)objective(x, W_base, rho), x0, [], [], [], [], lb, ub, ...
    @(x)nonlcon(x, W_base, E, sigma_allow, disp_limit, node_coords, members), options);

%% Results
disp('----------------------------------');
disp(['Optimal thickness t = ', num2str(x_opt(1)), ' m']);
disp(['Optimal height/width ratio r = ', num2str(x_opt(2))]);
disp(['Minimum mass = ', num2str(mass_opt), ' kg']);
disp('----------------------------------');



% %% Final Contour Plot
% figure;
% hold on;
% contourf(T, R, MASS, 20, 'LineColor', 'none');
% colorbar;
% xlabel('Thickness t (m)');
% ylabel('Height-to-width ratio r');
% title('Mass Contour Plot with Optimization Path');
% 
% % Overlay Constraint Boundaries
% contour(T, R, CONSTRAINT_VIOLATION, [1 1], 'k-', 'LineWidth', 2);
% 
% % Overlay Optimization Path
% if ~isempty(opt_path)
%     plot(opt_path(:,1), opt_path(:,2), 'wo-', 'MarkerFaceColor', 'w', 'LineWidth', 2);
% end
% 
% % Plot Final Optimized Point
% plot(x_opt(1), x_opt(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% text(x_opt(1) + 0.01, x_opt(2), 'Optimal Point', 'FontSize', 10, 'Color', 'r');
% 
% hold off;