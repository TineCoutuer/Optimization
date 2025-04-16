rho = 7850;        % kg/m^3
E = 210e9;         % Pa
sigma_allow = 250e6; % Pa
W_base = 0.3;      % m (constant width)
disp_limit = 0.02; % m (20 mm max displacement)
% Recommended Width: 0.05 m - 0.5m
% Recommended Height: 0.1m - 2.0 m
lb = [0.003, 1];       % Lower bounds [t_min, r_min]
ub = [0.05, 10];        % Upper bounds [t_max, r_max]
safety_fac = 1;     % Safety factor
F_ref = 1800000; %N (just to scale constraints)
L = 7; %delete

%% Nodes and Members
node_coords = [...
     0    0;
     7    0;
    14    0;
    21    0;
    3.5   7;
   10.5   7;
   17.5   7];

members = [
    1 2; 2 3; 3 4;       % bottom chord
    1 5; 5 2;            % left triangle
    2 6; 6 3;            % center triangle
    3 7; 7 4;            % right triangle
    5 6; 6 7             % top chord
];

