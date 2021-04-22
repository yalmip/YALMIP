function tests = test_global_ex2_1_9
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex2_1_9.gms
% Created 28-Jul-2007 18:43:44 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);
x5 = sdpvar(1);
x6 = sdpvar(1);
x7 = sdpvar(1);
x8 = sdpvar(1);
x9 = sdpvar(1);
x10 = sdpvar(1);
objvar = sdpvar(1);

% Define constraints 
F = ([]);
F=[F,-(x1*x2+x2*x3+x3*x4+x4*x5+x5*x6+x6*x7+x7*x8+x8*x9+x9*x10+x1*x3+x2*x4+x3*x5+x4*x6+x5*x7+x6*x8+x7*x9+x8*x10+x1*x9+x1*x10+x2*x10+x1*x5+x4*x7)-objvar==0];
F=[F,x1+x2+x3+x4+x5+x6+x7+x8+x9+x10==1];
F=[F,cut((x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)^2==1)];
F=[F,0<=x1];
F=[F,0<=x2];
F=[F,0<=x3];
F=[F,0<=x4];
F=[F,0<=x5];
F=[F,0<=x6];
F=[F,0<=x7];
F=[F,0<=x8];
F=[F,0<=x9];
F=[F,0<=x10];

% Solve problem
sol = optimize(F,objvar,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem == 0)
assert(abs(value(objvar)--0.375) <= 1e-2) 
