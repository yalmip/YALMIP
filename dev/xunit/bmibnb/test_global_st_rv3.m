function tests = test_global_st_rv3
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_rv3.gms
% Created 06-Aug-2007 09:17:27 using YALMIP R20070725

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
objective = -(-(-0.00055*sqr(x1)-0.0583*x1-0.0019*sqr(x2)-0.2318*x2-0.0002*sqr(x3)-0.0108*x3-0.00095*sqr(x4)-0.1634*x4-0.0046*sqr(x5)-0.138*x5-0.0035*sqr(x6)-0.357*x6-0.00315*sqr(x7)-0.1953*x7-0.00475*sqr(x8)-0.361*x8-0.0048*sqr(x9)-0.1824*x9-0.003*sqr(x10)-0.162*x10-0.00265*sqr(x11)-0.4346*x11-0.0017*sqr(x12)-0.1054*x12-0.0012*sqr(x13)-0.2376*x13-0.00295*sqr(x14)-0.0059*x14-0.00315*sqr(x15)-0.189*x15-0.0021*sqr(x16)-0.0252*x16-0.00225*sqr(x17)-0.099*x17-0.0034*sqr(x18)-0.3604*x18-0.001*sqr(x19)-0.022*x19-0.00305*sqr(x20)-0.3294*x20));

% Define constraints 
F = ([]);
F=[F,8*x1+5*x6+4*x7+6*x12+6*x13+9*x14+5*x19+x20<=220];
F=[F,3*x1+4*x2+3*x7+7*x8+4*x13+9*x14+3*x15+2*x20<=175];
F=[F,2*x2+x3+6*x8+8*x9+9*x14+9*x15+8*x16<=215];
F=[F,7*x3+x4+7*x9+9*x10+2*x15+4*x16+9*x17<=195];
F=[F,4*x4+4*x5+x10+3*x11+7*x16+2*x17+8*x18<=145];
F=[F,9*x5+5*x6+5*x11+7*x12+x17+4*x18+6*x19<=185];
F=[F,5*x1+5*x6+3*x7+8*x12+5*x13+9*x18+9*x19+x20<=225];
F=[F,x1+9*x2+9*x7+3*x8+9*x13+7*x14+4*x19+x20<=215];
F=[F,3*x1+6*x2+3*x3+4*x8+2*x9+6*x14+3*x15+8*x19+x20<=175];
F=[F,x2+2*x3+8*x4+4*x9+x10+9*x15+6*x16<=155];
F=[F,9*x3+3*x4+6*x5+x10+6*x11+9*x16+8*x17<=210];
F=[F,6*x4+3*x5+3*x6+6*x11+3*x12+8*x17+9*x18<=190];
F=[F,9*x5+8*x6+2*x7+7*x12+8*x13+4*x18+3*x19<=205];
F=[F,4*x1+6*x6+9*x7+x8+6*x13+9*x14+8*x19+6*x20<=245];
F=[F,7*x1+3*x2+7*x7+4*x8+2*x9+x14+3*x15+5*x20<=160];
F=[F,7*x2+9*x3+7*x8+9*x9+5*x10+2*x15+6*x16<=225];
F=[F,6*x3+9*x4+8*x9+4*x10+2*x11+6*x16+4*x17<=195];
F=[F,5*x4+5*x5+7*x10+8*x11+9*x12+8*x17+6*x18<=240];
F=[F,7*x5+5*x6+6*x11+2*x12+8*x13+6*x18+9*x19<=215];
F=[F,x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20<=200];
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
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--35.7607 ) <= 1e-1) 