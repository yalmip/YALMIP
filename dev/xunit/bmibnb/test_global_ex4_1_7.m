function test_global_ex4_1_7
% Model generated from ex4_1_7.gms
% Created 28-Jul-2007 18:54:22 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);

% Define objective function 
objective = -(-(power(x1,4)-3*power(x1,3)-1.5*sqr(x1)+10*x1));

% Define constraints 
F = ([]);
F=[F,-5<=x1<=5];

% Solve problem
sol = solvesdp(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));

mbg_asserttrue(sol.problem==0);
mbg_asserttolequal(double(objective), -7.5, 1e-2);