function tests = test_global_ex2_1_4
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex2_1_4.gms
% Created 28-Jul-2007 18:43:35 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);
x5 = sdpvar(1);
x6 = sdpvar(1);
objvar = sdpvar(1);

% Define constraints 
F = ([]);
F=[F,-(6.5*x1-0.5*x1*x1)+x2+2*x3+3*x4+2*x5+x6+objvar==0];
F=[F,x1+2*x2+8*x3+x4+3*x5+5*x6<=16];
F=[F,-8*x1-4*x2-2*x3+2*x4+4*x5-x6<=-1];
F=[F,2*x1+0.5*x2+0.2*x3-3*x4-x5-4*x6<=24];
F=[F,0.2*x1+2*x2+0.1*x3-4*x4+2*x5+2*x6<=12];
F=[F,-0.1*x1-0.5*x2+2*x3+5*x4-5*x5+3*x6<=3];
F=[F,0<=x1<=1];
F=[F,0<=x2];
F=[F,0<=x3];
F=[F,0<=x4<=1];
F=[F,0<=x5<=1];
F=[F,0<=x6<=2];

% Solve problem
sol = optimize(F,objvar,sdpsettings('solver','bmibnb','allownonconvex',1));

assert(sol.problem==0)
assert(abs(value(objvar)--11) <= 1e-2)