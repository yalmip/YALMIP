function tests = test_global_st_ph10
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_ph10.gms
% Created 06-Aug-2007 09:31:31 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);

% Define objective function 
objective = -(-(3*x1-1.5*sqr(x1)-3.5*sqr(x2)+7*x2)+0-(0));

% Define constraints 
F = ([]);
F=[F,-2*x1+x2<=1];
F=[F,x1+2*x2<=7];
F=[F,x1+x2<=5];
F=[F,x1-2*x2<=2];
F=[F,0<=x1];
F=[F,x2<=0];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--10.5) <= 1e-2)