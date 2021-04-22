function tests = test_global_harker
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from harker.gms
% Created 17-Mar-2008 16:01:38 using YALMIP R20070810

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
x18 = sdpvar(1);
x19 = sdpvar(1);
x20 = sdpvar(1);

% Define objective function 
objective = -(19*x15-0.1*sqr(x15)-0.5*sqr(x18)-x18-0.005*sqr(x16)+27*x16-0.4*sqr(x19)-2*x19-0.15*sqr(x17)+30*x17-0.3*sqr(x20)-1.5*x20-(0.166666666666667*power(x1,3)+x1+0.0666666666666667*power(x2,3)+2*x2+0.1*power(x3,3)+3*x3+0.133333333333333*power(x4,3)+x4+0.1*power(x5,3)+2*x5+0.0333333333333333*power(x6,3)+x6+0.0333333333333333*power(x7,3)+x7+0.166666666666667*power(x8,3)+3*x8+0.0666666666666667*power(x9,3)+2*x9+0.333333333333333*power(x10,3)+x10+0.0833333333333333*power(x11,3)+2*x11+0.0666666666666667*power(x12,3)+2*x12+0.3*power(x13,3)+x13+0.266666666666667*power(x14,3)+3*x14))-0-(0);

% Define constraints 
F = ([]);
F=[F,x15+x16+x17-x18-x19-x20==0];
F=[F,-x1-x2+x5+x8-x15+x18==0];
F=[F,-x3+x11-x16+x19==0];
F=[F,-x4+x12-x17+x20==0];
F=[F,x1-x5-x6-x7+x9+x13==0];
F=[F,x2+x6-x8-x9-x10+x14==0];
F=[F,x3+x4+x7+x10-x11-x12-x13-x14==0];
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
F=[F,0<=x11];
F=[F,0<=x12];
F=[F,0<=x13];
F=[F,0<=x14];
F=[F,0<=x15];
F=[F,0<=x16];
F=[F,0<=x17];
F=[F,0<=x18];
F=[F,0<=x19];
F=[F,0<=x20];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1))
assert(sol.problem==0)
assert(abs(value(objective)--986.5135) <= 1)