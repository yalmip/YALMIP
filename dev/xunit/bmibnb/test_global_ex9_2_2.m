function tests = test_global_ex9_2_2
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex9_2_2.gms
% Created 28-Jul-2007 17:58:38 using YALMIP R20070725

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
x10 = sdpvar(1);
x11 = sdpvar(1);

% Define constraints 
objvar = x2*x2+(x3-10)*(x3-10);
F = ([]);
F=[F,x2<=15];
F=[F,-x2+x3<=0];
F=[F,-x2<=0];
F=[F,x2+x3+x4==20];
F=[F,-x3+x5==0];
F=[F,x3+x6==20];
F=[F,x4*x8==0];
F=[F,x5*x9==0];
F=[F,x6*x10==0];
F=[F,x7*x11==0];
F=[F,2*x2+4*x3+x8-x9+x10==60];
F=[F,0<=x2];
F=[F,0<=x3];
F=[F,0<=x4<=20];
F=[F,0<=x5<=20];
F=[F,0<=x6<=20];
F=[F,0<=x7<=20];
F=[F,0<=x8<=20];
F=[F,0<=x9<=20];
F=[F,0<=x10<=20];
F=[F,0<=x11<=20];

% Solve problem
sol = optimize(F,objvar,sdpsettings('solver','bmibnb','allownonconvex',1,'quadprog.Algorithm','interior-point-convex'));
assert(sol.problem==0)
assert(abs(value(objvar)-100) <= 1)