function tests = test_global_ex9_2_4
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex9_2_4.gms
% Created 28-Jul-2007 18:02:12 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
objvar = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);
x5 = sdpvar(1);
x6 = sdpvar(1);
x7 = sdpvar(1);
x8 = sdpvar(1);
x9 = sdpvar(1);

% Define constraints 
objvar = (0.5*x4-1)*(x4-2)+(0.5*x5-1)*(x5-2);
F = ([]);
F=[F,-x3+x4+x5==0];
F=[F,-x4+x6==0];
F=[F,-x5+x7==0];
F=[F,x6*x8==0];
F=[F,x7*x9==0];
F=[F,x2+x4-x8==0];
F=[F,x2-x9==-1];
F=[F,0<=x3];
F=[F,0<=x4];
F=[F,0<=x5];
F=[F,0<=x6<=200];
F=[F,0<=x7<=200];
F=[F,0<=x8<=200];
F=[F,0<=x9<=200];

% Solve problem
sol = optimize(F,objvar,sdpsettings('solver','bmibnb','allownonconvex',1,'quadprog.Algorithm','interior-point-convex'));
assert(sol.problem==0)
assert(abs(value(objvar)-0.5) <= 1e-3)