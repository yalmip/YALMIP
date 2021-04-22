function tests = test_global_st_ph20
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_ph20.gms
% Created 06-Aug-2007 09:30:28 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);

% Define objective function 
objective = -(-(15*x1-sqr(x1)-sqr(x2)-2*x2)-x3+0-(0));

% Define constraints 
F = ([]);
F=[F,-4*x1-3*x2+4*x3<=30];
F=[F,4*x1+9*x2-2*x3<=114];
F=[F,2*x2-x3<=8];
F=[F,2*x1+15*x2-8*x3<=64];
F=[F,x2<=14];
F=[F,-4*x1+3*x2-2*x3<=-18];
F=[F,4*x1-9*x2+4*x3<=-6];
F=[F,-6*x1+5*x2-4*x3<=-40];
F=[F,4*x1-9*x2-3*x3<=-132];
F=[F,0<=x1];
F=[F,0<=x2];
F=[F,0<=x3];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--158) <= 1e-2)