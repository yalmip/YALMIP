function tests = test_global_ex4_1_1
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex4_1_1.gms
% Created 28-Jul-2007 18:50:00 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);

% Define objective function 
objective = -(-(power(x1,6)-2.08*power(x1,5)+0.4875*power(x1,4)+7.1*power(x1,3)-3.95*sqr(x1)-x1)-(0.1));

% Define constraints 
F = ([]);
F=[F,-2<=x1<=11];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));

assert(sol.problem==0)
assert(abs(value(objective)- -7.48731236490236) <= 1e-2) 