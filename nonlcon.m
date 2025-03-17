function [c, ceq] = nonlcon(x, W_base, E, sigma_allow, disp_limit, node_coords, members)
    t = x(1);
    r = x(2);

    W = W_base;
    H = r * W_base;

    A = 2*W*t + (H-2*t)*t;

    % Moment of inertia
    I_flange = (1/12)*W*t^3;
    I_web = (1/12)*t*(H-2*t)^3;
    d_NA = H/2 - t/2;
    A_segment = W*t;
    I_tot = 2*(I_flange + A_segment*d_NA^2) + I_web;

    % Run structural analysis for current t and r
    [axial_forces, U_global, member_lengths] = structural_analysis(E, A, node_coords, members);

    % Stress constraint
    stress_violation = max(abs(axial_forces)/A) - sigma_allow;

    % Buckling constraint (check only compression members)
    K_factor = 1;
    F_cr = (pi^2 * E * I_tot) ./ (K_factor * member_lengths).^2;
    compression_forces = axial_forces(axial_forces < 0);
    member_lengths_compression = member_lengths(axial_forces < 0);
    F_cr_compression = (pi^2 * E * I_tot) ./ (K_factor * member_lengths_compression).^2;
    buckling_violation = max(abs(compression_forces) - F_cr_compression);

    % Displacement constraint (check vertical at node 6, P6)
    Uy_P6 = U_global((5-1)*2 + 2);
    displacement_violation = abs(Uy_P6) - disp_limit;

    % Inequality constraints (c <= 0 is satisfied)
    c = [stress_violation;
         buckling_violation;
         displacement_violation];

    ceq = []; % No equality constraints
end
