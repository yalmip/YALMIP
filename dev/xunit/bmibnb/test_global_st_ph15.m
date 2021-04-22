function tests = test_global_st_ph15
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_ph15.gms
% Created 06-Aug-2007 09:31:22 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);

% Define objective function 
objective = -(-(16*x1-4*sqr(x1)-0.5*sqr(x2)+x2-2.5*sqr(x3)+15*x3-sqr(x4)+8*x4)+0-(0));

% Define constraints 
F = ([]);
F=[F,x1-x2+3*x3-2*x4<=6];
F=[F,-x1+4*x2+x3-2*x4<=7];
F=[F,2*x1+x2+2*x3+x4<=29];
F=[F,x1-x2+x3+x4<=11];
F=[F,0<=x1];
F=[F,0<=x2];
F=[F,0<=x3];
F=[F,0<=x4];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--392.7037) <= 1)