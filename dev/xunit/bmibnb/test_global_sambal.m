function tests = test_global_sambal
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from sambal.gms
% Created 22-Aug-2007 09:42:02 using YALMIP R20070810

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
x11 = sdpvar(1);
x12 = sdpvar(1);
x13 = sdpvar(1);
x14 = sdpvar(1);
x15 = sdpvar(1);
x16 = sdpvar(1);
x17 = sdpvar(1);

% Define objective function 
objective = -(-(0.0666666666666667*sqr(15-x1)+0.333333333333333*sqr(3-x2)+0.00769230769230769*sqr(130-x3)+0.0125*sqr(80-x4)+0.0666666666666667*sqr(15-x7)+0.00769230769230769*sqr(130-x8)+0.05*sqr(20-x9)+0.04*sqr(25-x10)+0.025*sqr(40-x11)+0.0181818181818182*sqr(55-x12)+0.00454545454545455*sqr(220-x13)+0.00526315789473684*sqr(190-x16)+0.00952380952380952*sqr(105-x17))+0-(0));

% Define constraints 
F = ([]);
F=[F,-x1-x2-x3-x4+x13==0];
F=[F,-x5+x14==0];
F=[F,-x6+x15==0];
F=[F,-x7-x8-x9+x16==0];
F=[F,-x10-x11-x12+x17==0];
F=[F,-x5-x6+x13==0];
F=[F,-x1-x7-x10+x14==0];
F=[F,-x2-x8-x11+x15==0];
F=[F,-x3-x12+x16==0];
F=[F,-x4-x9+x17==0];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1,'quadprog.Algorithm','interior-point-convex'))
assert(sol.problem==0)
assert(abs(value(objective)-3.968220361763031) <= 1e-2) 