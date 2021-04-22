function tests = test_global_st_bsj2
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_bsj2.gms
% Created 17-Mar-2008 11:00:00 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);

% Define objective function 
objective = -(-(2*x1-sqr(x1)-sqr(x2)-sqr(x3)+2*x3)+0-(0));

% Define constraints 
F = ([]);
F=[F,x1+x2-x3<=1];
F=[F,-x1+x2-x3<=-1];
F=[F,12*x1+5*x2+12*x3<=34.8];
F=[F,12*x1+12*x2+7*x3<=29.1];
F=[F,-6*x1+x2+x3<=-4.1];
F=[F,0<=x1];
F=[F,0<=x2];
F=[F,0<=x3];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)-1) <= 1e-2) 