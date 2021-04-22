function tests = test_global_st_rv7
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_rv7.gms
% Created 06-Aug-2007 09:16:48 using YALMIP R20070725

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
x21 = sdpvar(1);
x22 = sdpvar(1);
x23 = sdpvar(1);
x24 = sdpvar(1);
x25 = sdpvar(1);
x26 = sdpvar(1);
x27 = sdpvar(1);
x28 = sdpvar(1);
x29 = sdpvar(1);
x30 = sdpvar(1);

% Define objective function 
objective = -(-(-0.00165*sqr(x1)-0.1914*x1-0.0004*sqr(x2)-0.0384*x2-0.00285*sqr(x3)-0.3876*x3-0.00155*sqr(x4)-0.1116*x4-0.0038*sqr(x5)-0.4636*x5-0.0044*sqr(x6)-0.044*x6-0.0046*sqr(x7)-0.3588*x7-0.00085*sqr(x8)-0.0272*x8-0.00165*sqr(x9)-0.231*x9-0.0025*sqr(x10)-0.27*x10-0.00385*sqr(x11)-0.308*x11-0.00355*sqr(x12)-0.3692*x12-0.0015*sqr(x13)-0.288*x13-0.0037*sqr(x14)-0.407*x14-0.00125*sqr(x15)-0.1175*x15-0.00095*sqr(x16)-0.1045*x16-0.0048*sqr(x17)-0.1632*x17-0.0015*sqr(x18)-0.135*x18-0.0048*sqr(x19)-0.0864*x19-0.0007*sqr(x20)-0.1176*x20-0.0043*sqr(x21)-0.645*x21-0.0045*sqr(x22)-0.882*x22-0.00245*sqr(x23)-0.3283*x23-0.0004*sqr(x24)-0.0648*x24-0.0048*sqr(x25)-0.0864*x25-0.00485*sqr(x26)-0.4753*x26-0.00025*sqr(x27)-0.046*x27-0.00435*sqr(x28)-0.7917*x28-0.00365*sqr(x29)-0.7008*x29-0.0002*sqr(x30)-0.0384*x30)+0-0);

% Define constraints 
F = ([]);
F=[F,4*x1+7*x6+4*x7+8*x12+x13+3*x14+8*x19+6*x20+x25+8*x26<=425];
F=[F,7*x1+3*x2+7*x7+9*x8+9*x13+2*x14+6*x15+5*x20+7*x21+5*x26+8*x27<=450];
F=[F,7*x2+9*x3+8*x8+4*x9+3*x14+6*x15+4*x16+6*x21+5*x22+3*x27+2*x28<=380];
F=[F,6*x3+9*x4+7*x9+8*x10+8*x15+8*x16+6*x17+5*x22+3*x23+2*x28+x29<=415];
F=[F,5*x4+5*x5+6*x10+2*x11+9*x16+6*x17+9*x18+9*x23+3*x24+3*x29+4*x30<=360];
F=[F,7*x5+5*x6+6*x11+6*x12+8*x17+5*x18+x19+9*x24+6*x25+4*x30<=365];
F=[F,4*x1+5*x6+4*x7+4*x12+9*x13+6*x18+2*x19+2*x20+2*x25+x26<=300];
F=[F,2*x1+x2+3*x7+7*x8+9*x13+9*x14+x19+4*x20+6*x21+5*x26+5*x27<=370];
F=[F,9*x1+7*x2+x3+6*x8+8*x9+2*x14+4*x15+x20+4*x21+7*x22+2*x27+4*x28<=370];
F=[F,3*x2+4*x3+4*x4+7*x9+9*x10+7*x15+2*x16+3*x21+2*x22+2*x23+x28+8*x29<=320];
F=[F,8*x3+9*x4+5*x5+x10+3*x11+x16+4*x17+7*x22+6*x23+4*x24+2*x29+6*x30<=330];
F=[F,6*x4+5*x5+3*x6+5*x11+7*x12+9*x17+9*x18+4*x23+x24+6*x25+2*x30<=325];
F=[F,3*x5+9*x6+3*x7+8*x12+5*x13+4*x18+x19+3*x24+6*x25+5*x26<=285];
F=[F,6*x1+2*x6+4*x7+2*x8+9*x13+7*x14+8*x19+2*x20+8*x25+8*x26+6*x27<=425];
F=[F,x1+2*x2+x7+4*x8+x9+6*x14+3*x15+7*x20+6*x21+5*x26+7*x27+3*x28<=335];
F=[F,9*x2+3*x3+2*x8+x9+6*x10+9*x15+6*x16+7*x21+6*x22+7*x27+5*x28+5*x29<=415];
F=[F,6*x3+3*x4+5*x9+6*x10+3*x11+9*x16+8*x17+7*x22+4*x23+7*x28+x29+6*x30<=390];
F=[F,9*x4+8*x5+2*x10+7*x11+8*x12+8*x17+9*x18+2*x23+x24+7*x29+3*x30<=410];
F=[F,6*x5+9*x6+9*x11+6*x12+9*x13+4*x18+3*x19+3*x24+x25+9*x30<=370];
F=[F,x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22+x23+x24+x25+x26+x27+x28+x29+x30<=400];
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
F=[F,0<=x21];
F=[F,0<=x22];
F=[F,0<=x23];
F=[F,0<=x24];
F=[F,0<=x25];
F=[F,0<=x26];
F=[F,0<=x27];
F=[F,0<=x28];
F=[F,0<=x29];
F=[F,0<=x30];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--138.1875) <= 1)