function tests = test_global_st_e12
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_e12.gms
% Created 21-Aug-2007 18:42:14 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);

% Define objective function 
objective = -(-(x1^0.6+x2^0.6-6*x1)+4*x3-3*x4+0-(0));

% Define constraints 
F = ([]);
F=[F,-3*x1+x2-3*x3==0];
F=[F,x1+2*x3<=4];
F=[F,x2+2*x4<=4];
F=[F,0<=x1<=3];
F=[F,0<=x2<=4];
F=[F,0<=x3<=2];
F=[F,0<=x4<=1];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--4.51420165136193) <= 1e-2) 