function tests = test_global_st_rv2
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_rv2.gms
% Created 06-Aug-2007 09:17:59 using YALMIP R20070725

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
objective = -(-(-0.00015*sqr(x1)-0.0051*x1-0.00245*sqr(x2)-0.2205*x2-0.00095*sqr(x3)-0.0171*x3-0.0038*sqr(x4)-0.6384*x4-0.0029*sqr(x5)-0.435*x5-0.0024*sqr(x6)-0.4704*x6-0.0034*sqr(x7)-0.4556*x7-0.0018*sqr(x8)-0.2916*x8-0.00305*sqr(x9)-0.0549*x9-0.00025*sqr(x10)-0.0245*x10-0.00195*sqr(x11)-0.3588*x11-0.0008*sqr(x12)-0.1456*x12-0.0035*sqr(x13)-0.672*x13-0.0027*sqr(x14)-0.5184*x14-0.002*sqr(x15)-0.016*x15-0.0026*sqr(x16)-0.1404*x16-0.0048*sqr(x17)-0.2592*x17-0.00275*sqr(x18)-0.418*x18-0.00235*sqr(x19)-0.1081*x19-0.00275*sqr(x20)-0.264*x20));

% Define constraints 
F = ([]);
F=[F,6*x1+2*x2+4*x3+3*x5+4*x6+9*x7+5*x9+x10+9*x11+6*x12+7*x14+9*x15+2*x16+8*x18+2*x19+4*x20<=405];
F=[F,6*x1+5*x2+x3+8*x4+4*x6+3*x7+9*x8+6*x10+4*x11+7*x12+5*x13+2*x15+5*x16+8*x17+9*x19+8*x20<=450];
F=[F,8*x2+6*x3+2*x4+6*x5+4*x7+4*x8+6*x9+9*x11+4*x12+6*x13+9*x14+9*x16+9*x17+3*x18+x20<=430];
F=[F,8*x1+7*x3+3*x4+2*x5+x6+7*x8+4*x9+7*x10+3*x12+4*x13+x14+6*x15+2*x17+8*x18+9*x19<=360];
F=[F,x1+5*x2+5*x4+5*x5+x6+3*x7+5*x9+7*x10+4*x11+6*x13+x14+3*x15+4*x16+3*x18+5*x19+5*x20<=315];
F=[F,x1+8*x2+7*x3+x5+6*x6+x7+6*x8+7*x10+3*x11+6*x12+4*x14+6*x15+x16+4*x17+x19+4*x20<=330];
F=[F,5*x2+8*x3+7*x4+3*x6+3*x7+8*x8+6*x9+6*x11+4*x12+3*x13+4*x15+2*x16+5*x17+2*x18+4*x20<=350];
F=[F,x1+3*x3+2*x4+7*x5+2*x7+x8+x9+7*x10+4*x12+3*x13+5*x14+3*x16+6*x17+3*x18+x19<=245];
F=[F,5*x1+5*x2+2*x4+x5+9*x6+7*x8+4*x9+8*x10+5*x11+2*x13+4*x14+4*x15+4*x17+8*x18+9*x19+x20<=390];
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
assert(abs(value(objective)--64.4807 ) <= 1e-1) 