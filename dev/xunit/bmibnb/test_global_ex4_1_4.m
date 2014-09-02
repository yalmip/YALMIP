function test_global_ex4_1_4
% Model generated from ex4_1_4.gms
% Created 28-Jul-2007 18:51:43 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);

% Define objective function 
objective = -(-(4*sqr(x1)-4*power(x1,3)+power(x1,4)));

% Define constraints 
F = ([]);
F=[F,-5<=x1<=5];

% Solve problem
sol = solvesdp(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));

mbg_asserttrue(sol.problem==0);
mbg_asserttolequal(double(objective), 0, 1e-2);