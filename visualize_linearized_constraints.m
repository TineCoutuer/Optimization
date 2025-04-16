function visualize_linearized_constraints(x0)
    % VISUALIZE_LINEARIZED_CONSTRAINTS(x0)
    % Plots real and linearized constraints & objective around point x0 = [t, r]

    clf; hold off; format long;

    % Load problem parameters
    params;

    % Grid of design variables
    t = 0.0005:0.001:0.10;
    r = 0.1:0.1:10.0;
    [T, R] = meshgrid(t, r);              % Grid: T(i,j), R(i,j)
    n_r = length(r);
    n_t = length(t);

    % Prepare grid to evaluate objective and constraints
    points = [T(:), R(:)];                % Flattened (N × 2) matrix

    % Allocate outputs
    fobj     = zeros(n_r, n_t);
    flingrid = zeros(n_r, n_t);
    g        = zeros(n_r, n_t, 8);
    glin     = zeros(n_r, n_t, 8);

    % Objective and constraint evaluation
    g0 = nonlcon(x0, W_base, E, L, sigma_allow, disp_limit, ...
                 F_ref, node_coords, members, safety_fac);  % [9×1]
    dG = dgw(x0);  % [9×2]

    for idx = 1:size(points, 1)
        x = points(idx, :);

        % Real objective and constraint
        fobj(idx) = objective(x, W_base, rho);
        gval = nonlcon(x, W_base, E, L, sigma_allow, disp_limit, ...
                       F_ref, node_coords, members, safety_fac);
        for k = 1:8
            g(idx + (k-1)*n_r*n_t) = gval(k);
        end

        % Linearized objective (not used in plotting anymore)
        flingrid(idx) = flin_func(x, x0);
    end

    % True linearized constraint surfaces
    for k = 1:8
        glin_vals = g0(k) + (points - x0) * dG(k, :)';  % [N×1]
        glin(:, :, k) = reshape(glin_vals, n_r, n_t);   % [n_r × n_t]
    end

    % Reshape fobj for contouring
    fobj = reshape(fobj, n_r, n_t);
    flingrid = reshape(flingrid, n_r, n_t);
    for k = 1:8
        g(:, :, k) = reshape(g(:, :, k), n_r, n_t);
    end

    % Plot real objective contours
    contour(t, r, fobj, 'LineColor', [0.8 0.8 0.8]); hold on; grid on;
    contour(t, r, flingrid,'--', 'LineWidth', 1); % linearized objective
    % Constraint colors
    constraint_colors = {'r', 'b', 'm', 'g', 'c', [0.5 0 0.5], [1 0.5 0], 'y'};

    % Plot nonlinear and linearized constraints
    for k = 1:8
        contour(t, r, g(:, :, k), [0 0], 'Color', constraint_colors{k}, 'LineWidth', 0.7);
        contour(t, r, glin(:, :, k), [0 0], '--', 'Color', constraint_colors{k}, 'LineWidth', 2);
    end

    % Mark linearization point
    plot(x0(1), x0(2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

    % Labels and legend
    xlabel('Thickness t (m)');
    ylabel('Height/Width ratio r');
    title('Real (-) and Linearized (--) Constraints & Objective');

    legend('Mass contours','Mass contours lin', ...
        'g1', 'g1 lin', ...
        'g2', 'g2 lin', ...
        'g3', 'g3 lin', ...
        'g4', 'g4 lin', ...
        'g5', 'g5 lin', ...
        'g6', 'g6 lin', ...
        'g7', 'g7 lin', ...
        'g8', 'g8 lin', ...
        'Linearization point');
end


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
end

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
end