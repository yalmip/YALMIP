function tests = test_global_ex5_2_2_case3
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex5_2_2_case3.gms
% Created 28-Jul-2007 19:03:11 using YALMIP R20070725

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
objvar = sdpvar(1);

% Define constraints 
objvar = -9*x1-15*x2+6*x3+13*x4+10*x5+10*x6;
F = ([]);
%F=[F,-9*x1-15*x2+6*x3+13*x4+10*x5+10*x6-objvar==0];
F=[F,-x3-x4+x8+x9==0];
F=[F,x1-x5-x8==0];
F=[F,x2-x6-x9==0];
F=[F,x7*x8-2.5*x1+2*x5<=0];
F=[F,x7*x9-1.5*x2+2*x6<=0];
F=[F,x7*x8+x7*x9-3*x3-x4==0];
F=[F,0<=x1<=100];
F=[F,0<=x2<=200];
F=[F,0<=x3<=500];
F=[F,0<=x4<=500];
F=[F,0<=x5<=500];
F=[F,0<=x6<=500];
F=[F,0<=x7<=500];
F=[F,0<=x8<=500];
F=[F,0<=x9<=500];

% Solve problem
sol = optimize(F,objvar,sdpsettings('solver','bmibnb','allownonconvex',1));

assert(sol.problem==0)
assert(abs(value(objvar)- -750) <= 1e-1)