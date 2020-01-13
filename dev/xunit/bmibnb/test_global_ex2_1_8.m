function tests = test_global_ex2_1_8
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex2_1_8.gms
% Created 28-Jul-2007 18:43:42 using YALMIP R20070725

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
objvar = sdpvar(1);

% Define constraints 
F = ([]);
F=[F,-(300*x1-7*x1*x1-4*x2*x2+270*x2-6*x3*x3+460*x3-8*x4*x4+800*x4-12*x5*x5+740*x5-9*x6*x6+600*x6-14*x7*x7+540*x7-7*x8*x8+380*x8-13*x9*x9+300*x9-12*x10*x10+490*x10-8*x11*x11+380*x11-4*x12*x12+760*x12-7*x13*x13+430*x13-9*x14*x14+250*x14-16*x15*x15+390*x15-8*x16*x16+600*x16-4*x17*x17+210*x17-10*x18*x18+830*x18-21*x19*x19+470*x19-13*x20*x20+680*x20-17*x21*x21+360*x21-9*x22*x22+290*x22-8*x23*x23+400*x23-4*x24*x24+310*x24)+objvar==0];
F=[F,x1+x2+x3+x4==8];
F=[F,x5+x6+x7+x8==24];
F=[F,x9+x10+x11+x12==20];
F=[F,x13+x14+x15+x16==24];
F=[F,x17+x18+x19+x20==16];
F=[F,x21+x22+x23+x24==12];
F=[F,x1+x5+x9+x13+x17+x21==29];
F=[F,x2+x6+x10+x14+x18+x22==41];
F=[F,x3+x7+x11+x15+x19+x23==13];
F=[F,x4+x8+x12+x16+x20+x24==21];
F=[F,0<=x1<=100];
F=[F,0<=x2<=100];
F=[F,0<=x3<=100];
F=[F,0<=x4<=100];
F=[F,0<=x5<=100];
F=[F,0<=x6<=100];
F=[F,0<=x7<=100];
F=[F,0<=x8<=100];
F=[F,0<=x9<=100];
F=[F,0<=x10<=100];
F=[F,0<=x11<=100];
F=[F,0<=x12<=100];
F=[F,0<=x13<=100];
F=[F,0<=x14<=100];
F=[F,0<=x15<=100];
F=[F,0<=x16<=100];
F=[F,0<=x17<=100];
F=[F,0<=x18<=100];
F=[F,0<=x19<=100];
F=[F,0<=x20<=100];
F=[F,0<=x21<=100];
F=[F,0<=x22<=100];
F=[F,0<=x23<=100];
F=[F,0<=x24<=100];

% Solve problem
sol = optimize(F,objvar,sdpsettings('solver','bmibnb','allownonconvex',1));

assert(sol.problem==0)
assert(abs(value(objvar)-1.563899998838327e+004) <= 1e-2) 