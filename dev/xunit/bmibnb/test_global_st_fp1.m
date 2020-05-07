function tests = test_global_st_fp1
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_fp1.gms
% Created 24-Jul-2007 14:16:00 using YALMIP R20070523

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);
x5 = sdpvar(1);

% Define objective function 
objective = -(-(42*x1-50*sqr(x1)-50*sqr(x2)+44*x2-50*sqr(x3)+45*x3-50*sqr(x4)+47*x4-50*sqr(x5)+47.5*x5));

% Define constraints 
F = ([]);
F=[F,20*x1+12*x2+11*x3+7*x4+4*x5<=40];
F=[F,0<=x1<=1];
F=[F,0<=x2<=1];
F=[F,0<=x3<=1];
F=[F,0<=x4<=1];
F=[F,0<=x5<=1];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--17) <= 1e-2) 