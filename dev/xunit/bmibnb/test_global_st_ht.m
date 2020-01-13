function tests = test_global_st_ht
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_ht.gms
% Created 06-Aug-2007 09:39:13 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);

% Define objective function 
objective = -(-(2.4*x1-sqr(x1)-sqr(x2)+1.2*x2)+0-(0));

% Define constraints 
F = ([]);
F=[F,-2*x1+x2<=1];
F=[F,x1+x2<=4];
F=[F,0.5*x1-x2<=1];
F=[F,0<=x1<=3];
F=[F,0<=x2<=2];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--1.6) <= 1e-2)