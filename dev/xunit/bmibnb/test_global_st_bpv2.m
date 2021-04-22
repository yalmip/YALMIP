function tests = test_global_st_bpv2
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_bpv2.gms
% Created 17-Mar-2008 10:59:04 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);

% Define objective function 
objective = -(-(x2*x3+x2+x3-x1*x3)+0-(0));

% Define constraints 
F = ([]);
F=[F,-4*x1-x2>=-12];
F=[F,3*x1-x2>=-1];
F=[F,4*x3-x4<=12];
F=[F,-x3-x4>=-8];
F=[F,4*x3-x4>=-1];
F=[F,0<=x1<=4];
F=[F,0<=x2<=4];
F=[F,0<=x3<=5];
F=[F,0<=x4<=5];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--8) <= 1e-2) 