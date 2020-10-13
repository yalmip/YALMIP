function tests = test_global_st_e37
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_e37.gms
% Created 06-Aug-2007 09:45:53 using YALMIP R20070725

% % Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x4 = sdpvar(1);
x5 = sdpvar(1);

% Define objective function 
objective = -(-(sqr(x4+x5-1.9837)+sqr(x4*exp(-x1)+x5*exp(-x2)-0.8393)+sqr(x4*exp(-2*x1)+x5*exp(-2*x2)-0.4305)+sqr(x4*exp(-3*x1)+x5*exp(-3*x2)-0.2441)+sqr(x4*exp(-4*x1)+x5*exp(-4*x2)-0.1248)+sqr(x4*exp(-5*x1)+x5*exp(-5*x2)-0.0981)+sqr(x4*exp(-6*x1)+x5*exp(-6*x2)-0.0549)+sqr(x4*exp(-7*x1)+x5*exp(-7*x2)-0.0174)+sqr(x4*exp(-8*x1)+x5*exp(-8*x2)-0.0249)+sqr(x4*exp(-9*x1)+x5*exp(-9*x2)-0.0154)+sqr(x4*exp(-10*x1)+x5*exp(-10*x2)-0.0127))+0-(0));

% Define constraints 
F = ([]);
F=[F,x1-x2<=0];
F=[F,0<=x1<=100];
F=[F,0<=x2<=100];
F=[F,1==x4];
F=[F,1==x5];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1,'bmibnb.maxiter',100));
assert(sol.problem == 0)
assert(abs(value(objective)-1.041720240470086e-03) <= 1e-2)