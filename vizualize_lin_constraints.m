% Code base from Exercise 7.1
% Comparison of nonlinear objective and constraints and the linearized
% constraints

% Initialization
clf, hold off
format long

% Design point, x0, where problem is linearized to be given outside springw7ex1.
x0 = [0.06, 8];

params;
t = 0.0005:0.001:0.15;
r = 0.1:0.1:10.0;

% Matrix of output values for combinations of design variables D and d: 
for j = 1:length(r)
    for i = 1:length(t)
        x = [t(i), r(j)];
 	     % Real objective function
        f = objective(x, W_base, rho);
        % Grid value of objective function:
        fobj(j,i) = f; 
        
        % Linearized objective function
        flin = flinw7(x,x0);
        % Grid value of objective function:
        flingrid(j,i) = flin; 
        
        
        % Real constraints:
        g = nonlcon(x, W_base, E, L, sigma_allow, disp_limit, F_ref, node_coords, members, safety_fac);
        % Grid values of constraints:
        g1(j,i) = g(1);    % stress constraint
        g2(j,i) = g(2);    % displacement tension constraint
        g3(j,i) = g(3);    % displacement compression constraint
        g4(j,i) = g(4);    % buckling x constraint
        g5(j,i) = g(5);    % buckling y constraint
        g6(j,i) = g(6);    % W/t const
        g7(j,i) = g(7);    % H/t const
        g8(j,i) = g(8);    % Wedge constraint
        
        % Linearized constraints:
        glin = glinw7(x,x0);
        % Grid values of linearized constraints:
        g1lin(j,i) = glin(1);    % Scaled length constraint
        g2lin(j,i) = glin(2);    % Scaled lowest force constraint
        g3lin(j,i) = glin(3);    % Scaled highest force constraint
        g4lin(j,i) = glin(4);    % Scaled shear stress constraint
        g5lin(j,i) = glin(5);    % Scaled frequency constraint

    end
end

% Contour plot of real spring problem
contour(D, d, fobj)
xlabel('Coil diameter D (m)'), ylabel('Wire diameter d (m)'), ...
   title('Real (-) and linearized (--) spring mass optimization problem')
hold on
contour(D, d, g1, [0.0 0.0],'k')
contour(D, d, g2, [0.0 0.0],'b')
contour(D, d, g3, [0.0 0.0],'c')
contour(D, d, g4, [0.0 0.0],'g')
contour(D, d, g5, [0.0 0.0],'r')

% Contour plot of linerized spring problem
contour(D, d, flingrid,'--')
contour(D, d, g1lin, [0.0 0.0],'k--')
contour(D, d, g2lin, [0.0 0.0],'b--')
contour(D, d, g3lin, [0.0 0.0],'c--')
contour(D, d, g4lin, [0.0 0.0],'g--')
contour(D, d, g5lin, [0.0 0.0],'r--')

grid

% Plot marker in initial design point:
plot(x0(1),x0(2),'o');

clear all

%end 