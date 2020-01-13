function tests = test_global_st_e01
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_e01.gms
% Created 21-Aug-2007 18:41:04 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);

% Define objective function 
objective = -(x1+x2+0-(0));

% Define constraints 
F = ([]);
F=[F,x1*x2<=4];
F=[F,0<=x1<=6];
F=[F,0<=x2<=4];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--6.66667) <= 1e-2) 