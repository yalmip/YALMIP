function tests = test_global_ex2_1_2
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex2_1_2.gms
% Created 28-Jul-2007 18:43:32 using YALMIP R20070725

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
F=[F,-(-0.5*(x1*x1+x2*x2+x3*x3+x4*x4+x5*x5)-10.5*x1-7.5*x2-3.5*x3-2.5*x4-1.5*x5)+10*x6+objvar==0];
F=[F,6*x1+3*x2+3*x3+2*x4+x5<=6.5];
F=[F,10*x1+10*x3+x6<=20];
F=[F,0<=x1<=1];
F=[F,0<=x2<=1];
F=[F,0<=x3<=1];
F=[F,0<=x4<=1];
F=[F,0<=x5<=1];
F=[F,0<=x6];

% Solve problem
sol = optimize(F,objvar,sdpsettings('solver','bmibnb','allownonconvex',1));

assert(sol.problem==0)
assert(abs(value(objvar)--213) <= 1e-2)