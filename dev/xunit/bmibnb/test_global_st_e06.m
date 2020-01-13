function tests = test_global_st_e06
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_e06.gms
% Created 21-Aug-2007 18:42:00 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);

% Define objective function 
objective = -(0-(0));

% Define constraints 
F = ([]);
F=[F,x3*x3-0.000169*x1*power(x2,3)==0];
F=[F,x1+x2+x3==50];
F=[F,-3*x1+x2==0];
F=[F,0<=x1<=12.5];
F=[F,0<=x2<=37.5];
F=[F,0<=x3<=50];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)-0) <= 1e-2) 