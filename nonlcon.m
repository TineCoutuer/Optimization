function [c, ceq] = nonlcon(x, W_base, E, L,sigma_allow, disp_limit, F_ref,node_coords, members, safety_fac)
    t = x(1);
    r = x(2);


    S = safety_fac;    

    W = W_base;
    H = r * W_base;

    A = 2*W*t + (H-2*t)*t;

    % Moment of inertia
    I_block     =   (1/12)* W*H^3;
    I_sub       =   (1/12)* (W-t)*(H-2*t)^3;
    I_x         =   I_block - I_sub;

    I_flange    =  (1/12)* t*      W^3;
    I_web       =  (1/12)* (H-2*t)*t^3;
    I_y         =  I_web + 2 * I_flange;

    % Run structural analysis for current t and r
    [axial_forces, U_global, member_lengths] = structural_analysis(E, A, node_coords, members);

    % Tension - Stress constraint
    %stress_violation = max(abs(axial_forces)/A) - sigma_allow;
    stress_tension = safety_fac * max(axial_forces) /(A *sigma_allow) - 1;
        
    % Tension - Displacement constraint
    displacement_tension = safety_fac * max(axial_forces)* L /(A *E*disp_limit) - 1;

    % Compression - Displacement constraint
    displacement_compression = - safety_fac * min(axial_forces)* L /(A *E*disp_limit) - 1;
    % Compression - Buckling constraint 
    K_factor = 1;
    compression_forces = axial_forces(axial_forces < 0);
    member_lengths_compression = member_lengths(axial_forces < 0);
    % (x-axis)

    F_cr_compression_x = (pi^2 * E * I_x) ./ (K_factor.^2 * member_lengths_compression.^2);
    buckling_x = 1/F_ref*(F_cr_compression_x - S *min(compression_forces));
    % (y-axis)
    F_cr_compression_y = (pi^2 * E * I_y) ./ (K_factor.^2 * member_lengths_compression.^2);
    buckling_y = 1/F_ref*(F_cr_compression_y - S *min(compression_forces));
    %buckling_violation = max(abs(compression_forces) - F_cr_compression);
    
    % Displacement constraint (check vertical at node 6, P6)
    % Uy_P6 = U_global((5-1)*2 + 2);
    % displacement_violation = abs(Uy_P6) - disp_limit;

    % Inequality constraints (c <= 0 is satisfied)
    % c = [stress_violation;
    %      buckling_violation;
    %      displacement_violation];

    c = [stress_tension;
     displacement_tension;
     displacement_compression;
     buckling_x;
     buckling_y];

    ceq = []; % No equality constraints
end
