function tests = test_global_ex8_1_3
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex8_1_3.gms
% Created 02-Aug-2007 11:04:35 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);

% Define objective function 
objective = -(-(1+sqr(1+x1+x2)*(19+3*sqr(x1)-14*x1+6*x1*x2-14*x2+3*sqr(x2)))*(30+sqr(2*x1-3*x2)*(18+12*sqr(x1)-32*x1-36*x1*x2+48*x2+27*sqr(x2))));

% Define constraints 
F = ([]);
% Solve problem
x = recover(objective);
sol = optimize(F+[-100<=x<=100],objective,sdpsettings('solver','bmibnb','allownonconvex',1))
assert(sol.problem==0 | sol.problem == 3)
assert(abs(value(objective)-3) <= 1e-2) 