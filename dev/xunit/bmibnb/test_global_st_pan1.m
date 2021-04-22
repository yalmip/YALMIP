function tests = test_global_st_pan1
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_pan1.gms
% Created 06-Aug-2007 09:34:07 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);

% Define objective function 
objective = -(-(1.25*x1-2.5*sqr(x1)-5*sqr(x2)+2.5*x2-7.5*sqr(x3)+5*x3)+0-(0));

% Define constraints 
F = ([]);
F=[F,10*x1+0.2*x2-0.1*x3<=11];
F=[F,-0.3*x1+9*x2+0.2*x3<=8];
F=[F,-0.1*x1+0.4*x2+11*x3<=12];
F=[F,6*x1+8*x2+9*x3<=18];
F=[F,0<=x1];
F=[F,0<=x2];
F=[F,0<=x3];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--5.2837) <= 1e-2)