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
    dF = dfw(x0);
     
    % Linearized objective function:
    flin = fx0 + dF*(x - x0)';
        
    %end 


function dF = dfw(x);
    % Derivatives of objective function 
    % Input:
    % x  : design point "[t, r]" for which derivatives are computed.
    params;
    % Output:
    % dF  : [1x2] gradient (row)vector "[dfdx1 dfdx2]" of objective function.
    
    % Note: Constant parameter values are read within the function springobj2. 
    % Forward finite diffence gradients of objective function and constraints
    
    % Finite diffence step
    hx = 1.0e-8;
    
    % Gradient of objective function
    fx = objective(x, W_base, rho);
    fx1plush = objective([x(1)+hx, x(2)], W_base, rho);
    fx2plush = objective([x(1)+hx, x(2)], W_base, rho);
    dfdx1 = (fx1plush - fx)/hx;
    dfdx2 = (fx2plush - fx)/hx;
    dF = [dfdx1 dfdx2];
    
    % end 