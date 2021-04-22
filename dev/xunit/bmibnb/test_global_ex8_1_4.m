function tests = test_global_ex8_1_4
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex8_1_4.gms
% Created 02-Aug-2007 11:06:05 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);

% Define objective function 
objective = -(-(12*sqr(x1)-6.3*power(x1,4)+power(x1,6)-6*x1*x2+6*sqr(x2)));

% Define constraints 
F = ([]);
% Solve problem
F = ([]);
x = recover(objective);
sol = optimize(F+[-100<=x<=100],objective,sdpsettings('solver','bmibnb','allownonconvex',1))

assert((sol.problem == 3) | (sol.problem == 0))
assert(abs(value(objective)-0) <= 1e-2)