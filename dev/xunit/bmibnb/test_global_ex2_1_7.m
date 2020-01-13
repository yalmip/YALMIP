function tests = test_global_ex2_1_7
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex2_1_7.gms
% Created 28-Jul-2007 18:43:40 using YALMIP R20070725

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
objvar = sdpvar(1);

% Define constraints
F = ([]);
F=[F,0.5*(sqr(x1-2)+2*sqr(x2-2)+3*sqr(x3-2)+4*sqr(x4-2)+5*sqr(x5-2)+6*sqr(x6-2)+7*sqr(x7-2)+8*sqr(x8-2)+9*sqr(x9-2)+10*sqr(x10-2)+11*sqr(x11-2)+12*sqr(x12-2)+13*sqr(x13-2)+14*sqr(x14-2)+15*sqr(x15-2)+16*sqr(x16-2)+17*sqr(x17-2)+18*sqr(x18-2)+19*sqr(x19-2)+20*sqr(x20-2))+objvar==0];
F=[F,-3*x1+7*x2-5*x4+x5+x6+2*x8-x9-x10-9*x11+3*x12+5*x13+x16+7*x17-7*x18-4*x19-6*x20<=-5];
F=[F,7*x1-5*x3+x4+x5+2*x7-x8-x9-9*x10+3*x11+5*x12+x15+7*x16-7*x17-4*x18-6*x19-3*x20<=2];
F=[F,-5*x2+x3+x4+2*x6-x7-x8-9*x9+3*x10+5*x11+x14+7*x15-7*x16-4*x17-6*x18-3*x19+7*x20<=-1];
F=[F,-5*x1+x2+x3+2*x5-x6-x7-9*x8+3*x9+5*x10+x13+7*x14-7*x15-4*x16-6*x17-3*x18+7*x19<=-3];
F=[F,x1+x2+2*x4-x5-x6-9*x7+3*x8+5*x9+x12+7*x13-7*x14-4*x15-6*x16-3*x17+7*x18-5*x20<=5];
F=[F,x1+2*x3-x4-x5-9*x6+3*x7+5*x8+x11+7*x12-7*x13-4*x14-6*x15-3*x16+7*x17-5*x19+x20<=4];
F=[F,2*x2-x3-x4-9*x5+3*x6+5*x7+x10+7*x11-7*x12-4*x13-6*x14-3*x15+7*x16-5*x18+x19+x20<=-1];
F=[F,2*x1-x2-x3-9*x4+3*x5+5*x6+x9+7*x10-7*x11-4*x12-6*x13-3*x14+7*x15-5*x17+x18+x19<=0];
F=[F,-x1-x2-9*x3+3*x4+5*x5+x8+7*x9-7*x10-4*x11-6*x12-3*x13+7*x14-5*x16+x17+x18+2*x20<=9];
F=[F,x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20<=40];
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
sol = optimize(F,objvar,sdpsettings('solver','bmibnb','allownonconvex',1));

assert(sol.problem==0)
assert(abs(value(objvar)--4150.4) <= 40)
