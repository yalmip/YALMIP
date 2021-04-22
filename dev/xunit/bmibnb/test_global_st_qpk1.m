function tests = test_global_st_qpk1
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_qpk1.gms
% Created 28-Jul-2007 18:15:36 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);

% Define objective function 
objective = -(-(2*x1-2*x1*x1+2*x1*x2+3*x2-2*x2*x2));

% Define constraints 
F = ([]);
F=[F,-x1+x2<=1];
F=[F,x1-x2<=1];
F=[F,-x1+2*x2<=3];
F=[F,2*x1-x2<=3];
F=[F,0<=x1];
F=[F,0<=x2];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));

assert(sol.problem==0)
assert(abs(value(objective)--3) <= 1e-2) 