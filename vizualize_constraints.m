% Code base from Exercise 6.1 and 2.2
% Visualization of truss mass optimization problem,
% Initialization
clf, hold off, clear
% Load parameters
params;
W_base = 0.2;
% Combinations of design variables t and r 
t = [0.00005:0.01:0.5];
r = [0.1:0.1:20.0];

% Matrix of output values for combinations of design variables D and d: 
for j=1:1:length(r)
  for i=1:1:length(t)
    % Assignment of design variables:
    x(1) = t(i);
    x(2) = r(j);
 	 % Objective function
    f = objective(x, W_base, rho);
    % Grid value of objective function:
    fobj(j,i) = f; 
     
    % Scaled constraints:
    [g,~] = nonlcon(x, W_base, E,L, sigma_allow, disp_limit, F_ref,node_coords, members,safety_fac);
    % Grid values of constraints:
    g1(j,i) = g(1);    % stress constraint
    g2(j,i) = g(2);    % displacement tension constraint
    g3(j,i) = g(3);    % displacement compression constraint
    g4(j,i) = g(4);    % buckling x constraint
    g5(j,i) = g(5);    % buckling y constraint
    g6(j,i) = g(6);    % W/t const
    g7(j,i) = g(7);    % H/t const
    g8(j,i) = g(8);    % Wedge constraint
    g9(j,i) = g(9);    % Area constraint


  end
end

% Plotting the stress constraint (g1) as a colored surface
figure;
surf(t, r, g1);
shading interp; % Smooth shading
colorbar; % Add a colorbar

% Add a thick black contour at g1 = 0
hold on;
contour3(t, r, g1, [0 0], 'k', 'LineWidth', 2);
hold off;

xlabel('Material thickness t (m)');
ylabel('Heigth/Width ratio r');
zlabel('Stress Constraint (g1)');
title('Stress Constraint (g1) Visualization');

% Plotting the displacement tension constraint (g2) as a colored surface
figure;
surf(t, r, g2);
shading interp;
colorbar;

hold on;
contour3(t, r, g2, [0 0], 'k', 'LineWidth', 2);
hold off;

xlabel('Material thickness t (m)');
ylabel('Heigth/Width ratio r');
zlabel('Displacement Tension Constraint (g2)');
title('Displacement Tension Constraint (g2) Visualization');

% Plotting the displacement compression constraint (g3) as a colored surface
figure;
surf(t, r, g3);
shading interp;
colorbar;

hold on;
contour3(t, r, g3, [0 0], 'k', 'LineWidth', 2);
hold off;

xlabel('Material thickness t (m)');
ylabel('Heigth/Width ratio r');
zlabel('Displacement Compression Constraint (g3)');
title('Displacement Compression Constraint (g3) Visualization');

% Plotting the buckling x constraint (g4) as a colored surface
figure;
surf(t, r, g4);
shading interp;
colorbar;

hold on;
contour3(t, r, g4, [0 0], 'k', 'LineWidth', 2);
hold off;

xlabel('Material thickness t (m)');
ylabel('Heigth/Width ratio r');
zlabel('Buckling x Constraint (g4)');
title('Buckling x Constraint (g4) Visualization');

% Plotting the buckling y constraint (g5) as a colored surface
figure;
surf(t, r, g5);
shading interp;
colorbar;

hold on;
contour3(t, r, g5, [0 0], 'k', 'LineWidth', 2);
hold off;

xlabel('Material thickness t (m)');
ylabel('Heigth/Width ratio r');
zlabel('Buckling y Constraint (g5)');
title('Buckling y Constraint (g5) Visualization');

% Plotting the Width-Thickness Ratio constraint (g6) as a colored surface
figure;
surf(t, r, g6);
shading interp;
colorbar;

hold on;
contour3(t, r, g6, [0 0], 'k', 'LineWidth', 2);
hold off;

xlabel('Material thickness t (m)');
ylabel('Heigth/Width ratio r');
zlabel('Width-Thickness Ratio Constraint (g6)');
title('Width-Thickness Ratio Constraint (g6) Visualization');

% Plotting the Height-Thickness Ratio constraint (g7) as a colored surface
figure;
surf(t, r, g7);
shading interp;
colorbar;

hold on;
contour3(t, r, g7, [0 0], 'k', 'LineWidth', 2);
hold off;

xlabel('Material thickness t (m)');
ylabel('Heigth/Width ratio r');
zlabel('Height-Thickness Ratio Constraint (g7)');
title('Height-Thickness Ratio Constraint (g7) Visualization');

% Plotting the wedge height constraint (g8) as a colored surface
figure;
surf(t, r, g8);
shading interp;
colorbar;

hold on;
contour3(t, r, g8, [0 0], 'k', 'LineWidth', 2);
hold off;

xlabel('Material thickness t (m)');
ylabel('Heigth/Width ratio r');
zlabel('Wedge Height Constraint (g8)');
title('Wedge Height Constraint (g8) Visualization');

% Plotting the area constraint (g8) as a colored surface
figure;
surf(t, r, g9);
shading interp;
colorbar;

hold on;
contour3(t, r, g9, [0 0], 'k', 'LineWidth', 2);
hold off;

xlabel('Material thickness t (m)');
ylabel('Heigth/Width ratio r');
zlabel('Pos. Area Constraint (g9)');
title('Pos. Area Constraint (g9) Visualization');