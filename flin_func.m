function flin = flin_func(x, x0);
% Based on Exercise 7.1 
% Input:
%   x  : ([1x2] row) design variables (t and r)
%   x0 : ([1x2] row) design point where objective is linearized.
% Output:
%   flin  : [1x1] scalar of linearized objective function value

% Assignment of designvariables
t = x(1);
r = x(2);

% Constant parameter values
params;

% Analysis of valve spring in linearization point x0:
smass = objective(x, W_base, rho);
 
% Objective function in linearization point:
fx0 = smass;
 
% Gradient of objective function in linearization point:
dF = dfw7ex1(x0);
 
% Linearized objective function:
flin = fx0 + dF*(x - x0)';
    
%end 