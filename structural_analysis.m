function [axial_forces, U_global, member_lengths] = structural_analysis(E, A, node_coords, members)
    DOF = 2;
    num_nodes = size(node_coords, 1);
    total_DOF = DOF * num_nodes;
    num_members = size(members, 1);

    K_global = zeros(total_DOF, total_DOF);

    member_lengths = zeros(num_members, 1);

    % Assemble global stiffness matrix
    for i = 1:num_members
        n1 = members(i, 1);
        n2 = members(i, 2);

        x1 = node_coords(n1, 1);
        y1 = node_coords(n1, 2);
        x2 = node_coords(n2, 1);
        y2 = node_coords(n2, 2);

        L = sqrt((x2 - x1)^2 + (y2 - y1)^2);
        cx = (x2 - x1) / L;
        cy = (y2 - y1) / L;

        k_local = (E * A / L) * ...
            [ cx*cx  cx*cy -cx*cx -cx*cy;
              cx*cy  cy*cy -cx*cy -cy*cy;
             -cx*cx -cx*cy  cx*cx  cx*cy;
             -cx*cy -cy*cy  cx*cy  cy*cy];

        dof_map = [DOF*(n1-1)+1, DOF*(n1-1)+2, DOF*(n2-1)+1, DOF*(n2-1)+2];

        K_global(dof_map, dof_map) = K_global(dof_map, dof_map) + k_local;

        member_lengths(i) = L;
    end

    % External loads
    F_global = zeros(total_DOF, 1);
    F_global((4-1)*2 + 2) = -700e3;   % Node 5 vertical load
    F_global((5-1)*2 + 2) = -1400e3;  % Node 6 vertical load
    F_global((6-1)*2 + 2) = -700e3;   % Node 7 vertical load

    % Boundary conditions
    fixed_DOFs = [1, 2, (4-1)*2 + 2]; % P1 (x,y) and P4 (y)
    free_DOFs = setdiff(1:total_DOF, fixed_DOFs);

    K_ff = K_global(free_DOFs, free_DOFs);
    F_f = F_global(free_DOFs);

    U_f = K_ff \ F_f;

    U_global = zeros(total_DOF, 1);
    U_global(free_DOFs) = U_f;

    % Compute axial forces
    axial_forces = zeros(num_members, 1);
    for i = 1:num_members
        n1 = members(i, 1);
        n2 = members(i, 2);

        x1 = node_coords(n1, 1);
        y1 = node_coords(n1, 2);
        x2 = node_coords(n2, 1);
        y2 = node_coords(n2, 2);

        L = member_lengths(i);
        cx = (x2 - x1) / L;
        cy = (y2 - y1) / L;

        dof_map = [DOF*(n1-1)+1, DOF*(n1-1)+2, DOF*(n2-1)+1, DOF*(n2-1)+2];
        u_member = U_global(dof_map);

        axial_forces(i) = (E * A / L) * [-cx, -cy, cx, cy] * u_member;
    end
end
