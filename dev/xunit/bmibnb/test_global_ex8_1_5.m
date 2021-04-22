function tests = test_global_ex8_1_5
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex8_1_5.gms
% Created 02-Aug-2007 11:07:14 using YALMIP R20070725

% Setup a clean YALMIP environment
yalmip('clear')

% Define all variables
x1 = sdpvar(1);
x2 = sdpvar(1);

% Define objective function
objective = -(-(4*sqr(x1)-2.1*power(x1,4)+0.333333333333333*power(x1,6)+x1*x2-4*sqr(x2)+4*power(x2,4)));

% Define constraints
F = ([]);
% Solve problem
x = recover(objective);
sol = optimize(F+[-100<=x<=100],objective,sdpsettings('solver','bmibnb','allownonconvex',1))

assert(sol.problem == 0)
assert(abs(value(objective)--1.032) <= 1e-2)
