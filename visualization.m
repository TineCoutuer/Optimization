% Code base from Exercise 6.1 and 2.2
% Visualization of truss mass optimization problem,


% Initialization
clf, hold off, clear

% Load parameters
params;

% Combinations of design variables t and r 
t = [0.00005:0.01:0.35];
r = [0.1:0.1:10.0];

%[axial_forces, U_global, member_lengths] = structural_analysis(E, A, node_coords, members);

% t = [0.005:0.001:0.03];
% r = [1.0:0.05:10.0];

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

% Contour plot of scaled spring problem

contour(t, r, fobj,[0.0 200.0 800.0 1400.0 2e3 4e3 8e3 16e3 24e3 32e3 50e3 70e3 10e4 12e4 16e4 20e4 25e4 40e4 80e4 12e5],'ShowText','off') % mass contour

hold on, grid on;

% Plot constraint boundaries (solid lines) and infeasible regions (dashed)
contour(t, r, g1, [0.0 0.0], 'r', 'LineWidth', 2)            % stress tension constraint
contour(t, r, g1, [8.0 8.0],'r--') % Infeasible regiony

contour(t, r, g2, [0.0 0.0], 'b', 'LineWidth', 2)            % displacement tension constraint
contour(t, r, g2, [8.0 8.0],'b--')   % Infeasible region

contour(t, r, g3, [0.0 0.0], 'm', 'LineWidth', 2)            % displacement compression constraint
contour(t, r, g3, [8.0 8.0],'m--')   % Infeasible region

contour(t, r, g4, [0.0 0.0], 'g', 'LineWidth', 2)            % buckling x constraint
contour(t, r, g4, [0.4 0.4],'g--')   % Infeasible region
 
contour(t, r, g5, [0.0 0.0], 'c', 'LineWidth', 2)            % buckling y constraint
contour(t, r, g5, [0.2 0.2],'c--')   % Infeasible region

% g6 - Width/Thickness ratio
contour(t, r, g6, [0.0 0.0], 'LineWidth', 2,'Color', [0.5 0 0.5])
contour(t, r, g6, [0.1 0.1], '--', 'Color', [0.5 0 0.5])

% g7 - Height/Thickness ratio
contour(t, r, g7, [0.0 0.0], 'LineWidth', 2,'Color', [1 0.5 0])   % orange
contour(t, r, g7, [0.1 0.1], '--', 'Color', [1 0.5 0])

% g8 - Minimum Wedge Height
contour(t, r, g8, [0.0 0.0], 'k', 'LineWidth', 2)
contour(t, r, g8, [0.1 0.1], 'k--')

% g9 - Area constraint
contour(t, r, g9, [0.0 0.0], 'y','LineWidth', 2)   % purple
contour(t, r, g9, [0.03 0.03], 'y--')


% Labels and legend
xlabel('Material thickness t (m)'), ylabel('Heigth/Width ratio r'), ...
   title('Figure 1: Truss size optimization problem')
legend('smass contours', ...
    'stress', 'infeasible', ...
    'displacement tension', 'infeasible', ...
    'displacement compression', 'infeasible', ...
    'buckling x', 'infeasible', ...
    'buckling y', 'infeasible', ...
    'width/thickness', 'infeasible', ...
    'height/thickness', 'infeasible', ...
    'wedge height', 'infeasible', ...
    'positive area', 'infeasible');

hold off;



%end 