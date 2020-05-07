function tests = test_global_ex4_1_6
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex4_1_6.gms
% Created 28-Jul-2007 18:53:56 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);

% Define objective function 
objective = -(-(power(x1,6)-15*power(x1,4)+27*sqr(x1))-(250));

% Define constraints 
F = ([]);
F=[F,-5<=x1<=5];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));

assert(sol.problem==0)
assert(abs(value(objective)- 7) <= 1e-2) 