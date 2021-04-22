function tests = test_global_ramsey
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ramsey.gms
% Created 22-Aug-2007 09:42:37 using YALMIP R20070810

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
x31 = sdpvar(1);
x32 = sdpvar(1);
x33 = sdpvar(1);

% Define objective function
objective = -(0.95*log(x12)+0.9025*log(x13)+0.857375*log(x14)+0.81450625*log(x15)+0.7737809375*log(x16)+0.735091890625*log(x17)+0.69833729609375*log(x18)+0.663420431289062*log(x19)+0.630249409724609*log(x20)+0.598736939238379*log(x21)+11.3760018455292*log(x22))-0-(0);

% Define constraints
F = ([]);
F=[F,0.759835685651593*x1^0.25-x12-x23==0];
F=[F,0.77686866556676*x2^0.25-x13-x24==0];
F=[F,0.794283468039448*x3^0.25-x14-x25==0];
F=[F,0.812088652256959*x4^0.25-x15-x26==0];
F=[F,0.830292969275008*x5^0.25-x16-x27==0];
F=[F,0.848905366318769*x6^0.25-x17-x28==0];
F=[F,0.867934991180342*x7^0.25-x18-x29==0];
F=[F,0.88739119671479*x8^0.25-x19-x30==0];
F=[F,0.907283545436972*x9^0.25-x20-x31==0];
F=[F,0.92762181422141*x10^0.25-x21-x32==0];
F=[F,0.948415999107521*x11^0.25-x22-x33==0];
F=[F,-x1+x2-x23==0];
F=[F,-x2+x3-x24==0];
F=[F,-x3+x4-x25==0];
F=[F,-x4+x5-x26==0];
F=[F,-x5+x6-x27==0];
F=[F,-x6+x7-x28==0];
F=[F,-x7+x8-x29==0];
F=[F,-x8+x9-x30==0];
F=[F,-x9+x10-x31==0];
F=[F,-x10+x11-x32==0];
F=[F,0.03*x11-x33<=0];
F=[F,3==x1];
F=[F,3<=x2];
F=[F,3<=x3];
F=[F,3<=x4];
F=[F,3<=x5];
F=[F,3<=x6];
F=[F,3<=x7];
F=[F,3<=x8];
F=[F,3<=x9];
F=[F,3<=x10];
F=[F,3<=x11];
F=[F,0.95<=x12];
F=[F,0.95<=x13];
F=[F,0.95<=x14];
F=[F,0.95<=x15];
F=[F,0.95<=x16];
F=[F,0.95<=x17];
F=[F,0.95<=x18];
F=[F,0.95<=x19];
F=[F,0.95<=x20];
F=[F,0.95<=x21];
F=[F,0.95<=x22];
F=[F,0.05==x23];
F=[F,0.05<=x24<=0.0575];
F=[F,0.05<=x25<=0.066125];
F=[F,0.05<=x26<=0.076044];
F=[F,0.05<=x27<=0.08745];
F=[F,0.05<=x28<=0.10057];
F=[F,0.05<=x29<=0.11565];
F=[F,0.05<=x30<=0.133];
F=[F,0.05<=x31<=0.15295];
F=[F,0.05<=x32<=0.17589];
F=[F,0.05<=x33<=0.20228];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--2.487467848213407) <= 1e-2)
