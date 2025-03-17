function mass = objective(x, W_base, rho)
    t = x(1);
    r = x(2);

    W = W_base;
    H = r * W_base;

    A = 2*W*t + (H-2*t)*t;
       
    height = 6; %m
    L_A = 7;  % Horizontal members length
    L_B = sqrt((L_A/2)^2 + height^2);  % Diagonals depend on H

    mass = rho * A * (5*L_A + 6*L_B);
end
