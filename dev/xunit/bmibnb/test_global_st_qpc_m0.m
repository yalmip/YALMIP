function tests = test_global_st_qpc_m0
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_qpc-m0.gms
% Created 21-Aug-2007 18:36:53 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);

% Define objective function 
objective = -(-(2*x1-x1*x1-x2*x2+4*x2)+0-(0));

% Define constraints 
F = ([]);
F=[F,x1-4*x2>=-8];
F=[F,-3*x1+x2>=-9];
F=[F,0<=x1];
F=[F,0<=x2];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--5) <= 1e-2)