function plot_bridge_visualized()
    % Load system parameters
    params;

    figure; hold on; axis equal; grid on;
    title('Warren Truss Bridge with Forces and Axial Loads');
    xlabel('x [m]'); ylabel('y [m]');



    % Plot nodes
    plot(node_coords(:,1), node_coords(:,2), 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

    % Label nodes
    for i = 1:7
        offset = -0.7 * (i <= 4) + 0.7 * (i > 4);
        text(node_coords(i,1)-1, node_coords(i,2) + offset, ...
             sprintf('P%d', i), 'HorizontalAlignment', 'center');
    end

    % Structural analysis
    t_vis = 0.02; r_vis = 5; % just for visualization (not optimized)
    W = W_base; H = r_vis * W_base;
    A_vis = 2*W*t_vis + (H-2*t_vis)*t_vis;

    [axial_forces, ~, ~, R] = structural_analysis(E, A_vis, node_coords, members);

    % Bearing reaction forces at supports (P1: node 1, P4: node 4)
    % DOFs: P1 (x=1, y=2), P4 y=8
    % Note: Only y-direction at P4 (roller)
    support_info = {
        1, 1, 'x';  % P1 x
        1, 2, 'y';  % P1 y
        4, 2, 'y'   % P4 y
    };
    
    scale_R = 4e5;  % scale like loads

    for i = 1:size(support_info,1)
        n = support_info{i,1};
        dof = support_info{i,2};
        dir = support_info{i,3};
    
        x0 = node_coords(n,1);
        y0 = node_coords(n,2);
        global_dof = (n - 1) * 2 + dof;
        Fval = R(global_dof);
    
        dx = strcmp(dir, 'x') * Fval / scale_R;
        dy = strcmp(dir, 'y') * Fval / scale_R;
    
        quiver(x0, y0-dy, dx, dy, 0, 'g', 'LineWidth', 1.5, 'MaxHeadSize', 2);
        text(x0 + dx, y0 - dy - 0.5, sprintf('%.0f kN', Fval/1e3), ...
             'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'g');
    end

    % Plot members with color for tension/compression and annotate force
    for i = 1:size(members, 1)
        n1 = members(i, 1);
        n2 = members(i, 2);

        x_vals = node_coords([n1 n2], 1);
        y_vals = node_coords([n1 n2], 2);
        F = axial_forces(i);

        color = 'r'; % Tension
        if F < 0
            color = 'b'; % Compression
        end

        plot(x_vals, y_vals, '-', 'Color', color, 'LineWidth', 2);

        % Annotate axial force
        text(mean(x_vals)-1, mean(y_vals)-0.5, sprintf('%.0f kN', F/1e3), ...
            'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', color);
    end

    % Plot external loads as arrows (P5, P6, P7)
    top_nodes = [5, 6, 7];
    forces = [-700e3, -1400e3, -700e3];  % Node loads
    scale = 4e5;                         % scaling for quiver

    for i = 1:length(top_nodes)
        n = top_nodes(i);
        x0 = node_coords(n, 1);
        y0 = node_coords(n, 2);
        Fy = forces(i) / scale;

        quiver(x0, y0-Fy, 0, Fy, 0, 'k', 'MaxHeadSize', 2, 'LineWidth', 1.5);
        text(x0, y0 - Fy + 0.5, sprintf('%.0f kN', abs(forces(i))/1e3), ...
             'HorizontalAlignment', 'center', 'FontSize', 8);
    end

        % Dummy plots for legend
    h1 = plot(nan, nan, 'r-', 'LineWidth', 2); % Tension
    h2 = plot(nan, nan, 'b-', 'LineWidth', 2); % Compression
    h3 = plot(nan, nan, 'k-', 'LineWidth', 1.5); % External load
    h4 = plot(nan, nan, 'g-', 'LineWidth', 1.5); % Reaction Forces

    % Add legend
    legend([h1 h2 h3 h4], {'Tension', 'Compression', 'External Load','Reaction Forces'}, ...
           'Location', 'southoutside', 'Orientation', 'horizontal');

    grid off; xlim([-3, 24]); ylim([-6, 15])
end
