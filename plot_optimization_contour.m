function plot_optimization_contour(iter_points)
    clf, hold off, clearvars -except iter_points

    params;
    t = 0.0005:0.001:0.15;
    r = 0.1:0.1:10.0;

    for j = 1:length(r)
        for i = 1:length(t)
            x = [t(i), r(j)];
            fobj(j,i) = objective(x, W_base, rho);
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
            %g9(j,i) = g(9);    % Area constraint
        end
    end
    contour(t, r, fobj,[0.0 200.0 800.0 1400.0 2e3 4e3 8e3 16e3 24e3 32e3 50e3 70e3 10e4 12e4 16e4 20e4 25e4 40e4 80e4 12e5],'ShowText','off') % mass contour
    hold on; grid on;
    contour(t, r, g1, [0 0], 'r', 'LineWidth', 2); contour(t, r, g1, [0.1 0.1], 'r--');
    contour(t, r, g2, [0 0], 'b', 'LineWidth', 2); contour(t, r, g2, [0.1 0.1], 'b--');
    contour(t, r, g3, [0 0], 'm', 'LineWidth', 2); contour(t, r, g3, [0.1 0.1], 'm--');
    contour(t, r, g4, [0 0], 'g', 'LineWidth', 2); contour(t, r, g4, [0.1 0.1], 'g--');
    contour(t, r, g5, [0 0], 'c', 'LineWidth', 2); contour(t, r, g5, [0.1 0.1], 'c--');
    contour(t, r, g6, [0.0 0.0], 'LineWidth', 2,'Color', [0.5 0 0.5]); contour(t, r, g6, [0.1 0.1], '--', 'Color', [0.5 0 0.5])
    contour(t, r, g7, [0.0 0.0], 'LineWidth', 2,'Color', [1 0.5 0]); contour(t, r, g7, [0.1 0.1], '--', 'Color', [1 0.5 0])
    contour(t, r, g8, [0.0 0.0], 'y', 'LineWidth', 2); contour(t, r, g8, [0.1 0.1], 'y--')
    %contour(t, r, g9, [0.0 0.0], 'y','LineWidth', 2); contour(t, r, g9, [0.03 0.03], 'y--')

    xlabel('Material thickness t (m)');
    ylabel('Height/Width ratio r');
    title('Optimization Contour Plot with Iteration Path');


    % Optional: plot the path
    if exist('iter_points', 'var') && ~isempty(iter_points)
        plot(iter_points(:,1), iter_points(:,2), 'k.-', 'LineWidth', 2, 'MarkerSize', 12);
        plot(iter_points(end,1), iter_points(end,2), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
        text(iter_points(end,1), iter_points(end,2)+0.2, 'Final Optimum', 'Color', 'red');
    end


    legend('mass contours', ...
    'g1 stress', 'infeasible', ...
    'g2 displacement tension', 'infeasible', ...
    'g3 displacement compression', 'infeasible', ...
    'g4 buckling x', 'infeasible', ...
    'g5 buckling y', 'infeasible', ...
    'g6 width/thickness', 'infeasible', ...
    'g7 height/thickness', 'infeasible', ...
    'g8 wedge height', 'infeasible', ...
    'g9 positive area', 'infeasible',...
    'optimization iteration','Final Optimum');

        % Display constraint values at final point
    if exist('iter_points', 'var') && ~isempty(iter_points)
        x_final = iter_points(end,:);
        [c_final, ~] = nonlcon(x_final, W_base, E, L, sigma_allow, disp_limit, F_ref, node_coords, members, safety_fac);
        disp('--- Constraint values at final optimum ---');
        for i = 1:length(c_final)
            fprintf('g%d = %.4f\n', i, c_final(i));
        end
    end


    hold off;
end


