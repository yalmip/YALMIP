function tests = test_global_st_ph12
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_ph12.gms
% Created 06-Aug-2007 09:31:28 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);

% Define objective function 
objective = -(-(x1-0.5*sqr(x1)-sqr(x2)+2*x2-1.5*sqr(x3)+3*x3)+0-(0));

% Define constraints 
F = ([]);
F=[F,x1<=4];
F=[F,x2<=4];
F=[F,x3<=4];
F=[F,2*x1+3*x2+4*x3<=35];
F=[F,0<=x1];
F=[F,0<=x2];
F=[F,0<=x3];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--22.625) <= 1e-1)