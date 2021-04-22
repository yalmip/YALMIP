function tests = test_global_ex2_1_1
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex2_1_1.gms
% Created 28-Jul-2007 18:42:42 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);
x5 = sdpvar(1);
objvar = sdpvar(1);

% Define constraints 
F = ([]);
F=[F,-(42*x1-0.5*(100*x1*x1+100*x2*x2+100*x3*x3+100*x4*x4+100*x5*x5)+44*x2+45*x3+47*x4+47.5*x5)+objvar==0];
F=[F,20*x1+12*x2+11*x3+7*x4+4*x5<=40];
F=[F,0<=x1<=1];
F=[F,0<=x2<=1];
F=[F,0<=x3<=1];
F=[F,0<=x4<=1];
F=[F,0<=x5<=1];

% Solve problem
sol = optimize(F,objvar,sdpsettings('solver','bmibnb','allownonconvex',1));

assert(sol.problem==0)
assert(abs(value(objvar)--17) <= 1e-2) 