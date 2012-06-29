% Model generated from ex4_1_8.gms
% Created 28-Jul-2007 18:54:37 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
objvar = sdpvar(1);

% Define constraints 
objvar = -(-(sqr(x2)-7*x2)+12*x1);
F = set([]);
F=[F,-2*POWER(x1,4)-x2==-2];
F=[F,0<=x1<=2];
F=[F,0<=x2<=3];

% Solve problem
sol = solvesdp(F,objvar,sdpsettings('solver','bmibnb','allownonconvex',1));

mbg_assertfalse(sol.problem);
mbg_asserttolequal(double(objvar), -16.7389, 1e-2);