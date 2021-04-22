function tests = test_global_dispatch
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from dispatch.gms
% Created 28-Jul-2007 17:34:32 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);
objvar = sdpvar(1);

% Define constraints 
F = ([]);
F=[F,-(0.00533*sqr(x1)+11.669*x1+0.00889*sqr(x2)+10.333*x2+0.00741*sqr(x3)+10.833*x3)+objvar==653.1];
F=[F,-(0.01*(0.0676*x1*x1+0.00953*x1*x2-0.00507*x1*x3+0.00953*x2*x1+0.0521*x2*x2+0.00901*x2*x3-0.00507*x3*x1+0.00901*x3*x2+0.0294*x3*x3)-0.000766*x1-3.42e-5*x2+0.000189*x3)+x4==0.040357];
F=[F,x1+x2+x3-x4>=210];
F=[F,50<=x1<=200];
F=[F,37.5<=x2<=150];
F=[F,45<=x3<=180];

% Solve problem
sol = optimize(F,objvar,sdpsettings('solver','bmibnb','allownonconvex',1));

assert(sol.problem==0)
assert(abs(value(objvar)-3.1553e+003) <= 32)