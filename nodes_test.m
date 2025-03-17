% bridge optimization


clc; clear; close all;

%% Node Coordinates (x, y)
x = [0 7 14 21 3.5 10.5 17.5];
y = [0 0 0 0 6 6 6];

figure
hold on
axis equal
grid on

% Plot nodes
plot(x, y, 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 8)

% Add node labels
for i = 1:4
    text(x(i), y(i)-0.7, ['P' num2str(i)], 'HorizontalAlignment', 'center');
end

% Add node labels
for i = 5:length(x)
    text(x(i), y(i)+0.7, ['P' num2str(i)], 'HorizontalAlignment', 'center');
end

% Connect P1 with P2 and P5
plot([x(1), x(2)], [y(1), y(2)], 'b-')
plot([x(1), x(5)], [y(1), y(5)], 'b-')

% Connect P2 with P5 and P6 P3
plot([x(2), x(3)], [y(2), y(3)], 'b-')
plot([x(2), x(5)], [y(2), y(5)], 'b-')
plot([x(2), x(6)], [y(2), y(6)], 'b-')

% Connect P3 with P6 and P7 and P4
plot([x(3), x(4)], [y(3), y(4)], 'b-')
plot([x(3), x(6)], [y(3), y(6)], 'b-')
plot([x(3), x(7)], [y(3), y(7)], 'b-')

plot([x(4), x(7)], [y(4), y(7)], 'b-')
plot([x(5), x(6)], [y(5), y(6)], 'b-')
plot([x(6), x(7)], [y(6), y(7)], 'b-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters 
load = 200; %kN/m = combination of dead and traffic load is 200 kN/m.
height = 6; %m
width = 21; %m
length_hor = 7; %m

%% properties steel
rho = 7850; %kg/m^3 = DENSITY
E = 210000; %MPA = ELASTICITY / YOUNGS MODULUS
G = 81000; %MPA = SHEAR MODULUS


%% variables
% % length is based on distance nodes
% SEE OVERLEAF FOR WHAT LETTER IS WHAT DISTANCE
H = 50; % cm?
t = 5; % cm? 
W = 20; % cm
HW_ratio = H/W; % now it is = 2.5

%% Total mass of structure
% L1 = L2 = L3 = L4 = L5  = 7 m (straight bars)
L_A = 7; %m

% L6 = L7 = L8 = L9 = L10 = L11 = 7 m (slant bars)
L_B = sqrt((L_A/2)^2 + height^2);

% cross section of beam
A = 2*W*t + (H-2*t)*t;
% mass = sum(rho_i * A_i * L_i)
mass_i = rho * A;
mass_total = rho * (L_A*5 + L_B*6) * A*11;

% moment of intertia
I_1_3 = (1/12)*W*t^3; %  upper and lower segment = I_1 = I_3
I_2 = (1/12)*t*(H-2*t)^3; % middle segment

d_NA = H/2 - t/2; % vert distance from centroud of segment to neutral axis (same for upper and lower segment and zero for middle segment)
A_segment = W*t; % area of segment

I_tot = (I_1_3 + A_segment*d_NA^2) + (I_2) + (I_1_3 + A_segment*d_NA^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CONSTRAINTS

%% stress
% Normal stress in each member should be below allowable stress (F_I is
% force in member i)

% stress_i = F_i / A


%% buckling 
% for compression elements
% F_I should be less than cripling force!
% F_i <= (pi^2 * E*I) / (K*L_i)^2 

%% displacement?


%% constraints thickness and ratio?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STRUCTURAL ANALYSIS
% assuming P1 is pinned and P4 is roller

%% calculate load in each node
% node 1 and 6
F1_v = load * height ; %kN
F4_v = F1_v; %kN

F2_v = load * length_hor;
F3_v = F2_v;  %kN

%% Vertical load in endpoint
A_v = (F1_v + F2_v + F3_v + F4_v)/2;
B_v = A_v;



%% Calculate stress


%% Calculate displacement


