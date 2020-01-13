function tests = test_global_st_fp8
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_fp8.gms
% Created 06-Aug-2007 09:42:05 using YALMIP R20070725

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

% Define objective function 
objective = -(-(300*x1-7*sqr(x1)-4*sqr(x2)+270*x2-6*sqr(x3)+460*x3-8*sqr(x4)+800*x4-12*sqr(x5)+740*x5-9*sqr(x6)+600*x6-14*sqr(x7)+540*x7-7*sqr(x8)+380*x8-13*sqr(x9)+300*x9-12*sqr(x10)+490*x10-8*sqr(x11)+380*x11-4*sqr(x12)+760*x12-7*sqr(x13)+430*x13-9*sqr(x14)+250*x14-16*sqr(x15)+390*x15-8*sqr(x16)+600*x16-4*sqr(x17)+210*x17-10*sqr(x18)+830*x18-21*sqr(x19)+470*x19-13*sqr(x20)+680*x20-17*sqr(x21)+360*x21-9*sqr(x22)+290*x22-8*sqr(x23)+400*x23-4*sqr(x24)+310*x24)+0-(0));

% Define constraints 
F = ([]);
F=[F,-x1-x5-x9-x13-x17-x21<=-29];
F=[F,x1+x5+x9+x13+x17+x21<=29];
F=[F,-x2-x6-x10-x14-x18-x22<=-41];
F=[F,x2+x6+x10+x14+x18+x22<=41];
F=[F,-x3-x7-x11-x15-x19-x23<=-13];
F=[F,x3+x7+x11+x15+x19+x23<=13];
F=[F,-x4-x8-x12-x16-x20-x24<=-21];
F=[F,x4+x8+x12+x16+x20+x24<=21];
F=[F,-x1-x2-x3-x4<=-8];
F=[F,x1+x2+x3+x4<=8];
F=[F,-x5-x6-x7-x8<=-24];
F=[F,x5+x6+x7+x8<=24];
F=[F,-x9-x10-x11-x12<=-20];
F=[F,x9+x10+x11+x12<=20];
F=[F,-x13-x14-x15-x16<=-24];
F=[F,x13+x14+x15+x16<=24];
F=[F,-x17-x18-x19-x20<=-16];
F=[F,x17+x18+x19+x20<=16];
F=[F,-x21-x22-x23-x24<=-12];
F=[F,x21+x22+x23+x24<=12];
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

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)- 15639) <= 10)