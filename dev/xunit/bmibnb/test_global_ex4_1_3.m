function tests = test_global_ex4_1_3
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex4_1_3.gms
% Created 28-Jul-2007 18:50:55 using YALMIP R20070725

% Setup a clean YALMIP environment
yalmip('clear')

% Define all variables
x1 = sdpvar(1);

% Define objective function
objective = -(-(8.9248e-5*x1-0.0218343*sqr(x1)+0.998266*power(x1,3)-1.6995*power(x1,4)+0.2*power(x1,5)));

% Define constraints
F = ([]);
F=[F,0<=x1<=10];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1))
assert(sol.problem == 0)
assert(abs(value(objective)- -4.43672e2) <= 1e-2)
