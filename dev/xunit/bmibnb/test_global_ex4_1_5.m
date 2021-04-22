function tests = test_global_ex4_1_5
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex4_1_5.gms
% Created 28-Jul-2007 18:52:01 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);

% Define objective function 
objective = -(-(2*sqr(x1)-1.05*power(x1,4)+0.166666666666667*power(x1,6)-x1*x2+sqr(x2)));

% Define constraints 
F = ([]);
F=[F,-5<=x1];
F=[F,x2<=0];

% Solve problem
x = recover(F);
sol = optimize(F+[-100<=x<=100],objective,sdpsettings('solver','bmibnb','allownonconvex',1));

assert(sol.problem==0)
assert(abs(value(objective)- 0) <= 1e-2) 