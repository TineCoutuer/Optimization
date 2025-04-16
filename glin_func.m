function glin = glin(x, x0);
% Two variable valve spring problem - Exercise 7.1 
% Scaled linearized constraints. 
% Input:
%   x  : ([1x2] row) design variables (D and d)
%   x0 : ([1x2] row) design point where constraints are linearized.
% Output:
%   glin  : [1x5] row of linearized constraint values.

% Assignment of designvariables
D = x(1);
d = x(2);

% Constant parameter values
springparams2;

% Nominal valve area:
Av = Dv^2*pi/4;


% Analysis of valve spring in linearization point x0:
[svol,smass,bvol,matc,manc,Lmin,L2,k,F1,F2,Tau1,Tau2,freq1]=...
    springanalysis1(x0(1),x0(2),L0,L1,n,E,G,rho,Dv,h,p1,p2,nm,ncamfac,nne,matp,bldp);
 
% Scaled constraints linearization point:
 
% Scaled length constraint
g0(1) = Lmin/L2 - 1;
    
% Scaled lowest force constraint
F1min = Av * p1;
g0(2) = 1 - F1/F1min;
    
% Scaled highest force constraint
F2min = Av * p2;
g0(3) = 1 - F2/F2min;
    
% Scaled shear stress constraint
Tau12max = 600*10^6;
g0(4) = Tau2/Tau12max - 1;
    
% Scaled frequency constraint
freq1lb = ncamfac * nm/2;
g0(5) = 1 - freq1/freq1lb;
 
% Gradient of objective function in linearization point:
dG = dgw7ex1(x0);

% Linearized constraints:
glin = g0 + (dG*(x - x0)')';
    
%end 