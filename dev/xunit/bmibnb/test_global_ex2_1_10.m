function tests = test_global_ex2_1_10
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex2_1_10.gms
% Created 28-Jul-2007 18:43:46 using YALMIP R20070725

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
F=[F,-(0.5*(42*sqr(52+x11)+98*sqr(3+x12)+48*sqr(x13-81)+91*sqr(x14-30)+11*sqr(85+x15)+63*sqr(x16-68)+61*sqr(x17-27)+61*sqr(81+x18)+38*sqr(x19-97)+26*sqr(73+x20))-0.5*(63*sqr(19+x1)+15*sqr(27+x2)+44*sqr(23+x3)+91*sqr(53+x4)+45*sqr(42+x5)+50*sqr(x6-26)+89*sqr(33+x7)+58*sqr(23+x8)+86*sqr(x9-41)+82*sqr(x10-19)))+objvar==0];
F=[F,3*x1+5*x2+5*x3+6*x4+4*x5+4*x6+5*x7+6*x8+4*x9+4*x10+8*x11+4*x12+2*x13+x14+x15+x16+2*x17+x18+7*x19+3*x20<=380];
F=[F,5*x1+4*x2+5*x3+4*x4+x5+4*x6+4*x7+2*x8+5*x9+2*x10+3*x11+6*x12+x13+7*x14+7*x15+5*x16+8*x17+7*x18+2*x19+x20<=415];
F=[F,x1+5*x2+2*x3+4*x4+7*x5+3*x6+x7+5*x8+7*x9+6*x10+x11+7*x12+2*x13+4*x14+7*x15+5*x16+3*x17+4*x18+x19+2*x20<=385];
F=[F,3*x1+2*x2+6*x3+3*x4+2*x5+x6+6*x7+x8+7*x9+3*x10+7*x11+7*x12+8*x13+2*x14+3*x15+4*x16+5*x17+8*x18+x19+2*x20<=405];
F=[F,6*x1+6*x2+6*x3+4*x4+5*x5+2*x6+2*x7+4*x8+3*x9+2*x10+7*x11+5*x12+3*x13+6*x14+7*x15+5*x16+8*x17+4*x18+6*x19+3*x20<=470];
F=[F,5*x1+5*x2+2*x3+x4+3*x5+5*x6+5*x7+7*x8+4*x9+3*x10+4*x11+x12+7*x13+3*x14+8*x15+3*x16+x17+6*x18+2*x19+8*x20<=415];
F=[F,3*x1+6*x2+6*x3+3*x4+x5+6*x6+x7+6*x8+7*x9+x10+4*x11+3*x12+x13+4*x14+3*x15+6*x16+4*x17+6*x18+5*x19+4*x20<=400];
F=[F,x1+2*x2+x3+7*x4+8*x5+7*x6+6*x7+5*x8+8*x9+7*x10+2*x11+3*x12+5*x13+5*x14+4*x15+5*x16+4*x17+2*x18+2*x19+8*x20<=460];
F=[F,8*x1+5*x2+2*x3+5*x4+3*x5+8*x6+x7+3*x8+3*x9+5*x10+4*x11+5*x12+5*x13+6*x14+x15+7*x16+x17+2*x18+2*x19+4*x20<=400];
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
sol = optimize(F,objvar,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objvar)-49318) <= 150)