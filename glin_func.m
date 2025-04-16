function glin = glin_func(x, x0);
    % Two variable valve spring problem - Exercise 7.1 
    % Scaled linearized constraints. 
    % Input:
    %   x  : ([1x2] row) design variables (D and d)
    %   x0 : ([1x2] row) design point where constraints are linearized.
    % Output:
    %   glin  : [1x8] row of linearized constraint values.
    
    % Assignment of designvariables
    t = x(1);
    r = x(2);
    
    % Constant parameter values
    params;
    
    g0 = nonlcon(x, W_base, E, L, sigma_allow, disp_limit, F_ref, node_coords, members, safety_fac);
     
    % Scaled constraints linearization point:
     
    % Gradient of objective function in linearization point:
    dG = dgw(x0);
    
    % Linearized constraints:
    glin = g0(:)' + (dG * (x(:) - x0(:)))';
        
    %end 

function dG = dgw(x);
    % Derivatives of constraints
    params;
    % Input:
    % x  : design point "[t, r]" for which derivatives are computed.
    
    % Output:
    % dG  : [9X2] matrix with gradients of 5 constraints:
    %        "[dg1dx1  dg1dx2
    %          ....     ....
    %          dg9dx1  dg9dx2]" 
    
    % Note: Constant parameter values are read within the function
    
    % Forward finite diffence gradients of objective function and constraints
    
    % Finite diffence step
    hx = 1.0e-8;
    
    % Constraint gradients 
    gx = nonlcon(x, W_base, E, L, sigma_allow, disp_limit, F_ref, node_coords, members, safety_fac);
    gx1plush = nonlcon([x(1)+hx, x(2)], W_base, E, L, sigma_allow, disp_limit, F_ref, node_coords, members, safety_fac);
    gx2plush = nonlcon([x(1), x(2)+hx], W_base, E, L, sigma_allow, disp_limit, F_ref, node_coords, members, safety_fac);
    dgdx1 = (gx1plush - gx)./hx;
    dgdx2 = (gx2plush - gx)./hx;
    dG = [dgdx1, dgdx2];
    
    % end 
