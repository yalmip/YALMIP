function tests = test_global_hydro
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from hydro.gms
% Created 24-Jul-2007 13:32:20 using YALMIP R20070523

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

% Define objective function 
objective = -(-82.8*(0.0016*sqr(x1)+8*x1+0.0016*sqr(x2)+8*x2+0.0016*sqr(x3)+8*x3+0.0016*sqr(x4)+8*x4+0.0016*sqr(x5)+8*x5+0.0016*sqr(x6)+8*x6)-(248400));

% Define constraints 
F = ([]);
F=[F,x1+x7-x13>=1200];
F=[F,x2+x8-x14>=1500];
F=[F,x3+x9-x15>=1100];
F=[F,x4+x10-x16>=1800];
F=[F,x5+x11-x17>=950];
F=[F,x6+x12-x18>=1300];
F=[F,12*x19-x25+x26==24000];
F=[F,12*x20-x26+x27==24000];
F=[F,12*x21-x27+x28==24000];
F=[F,12*x22-x28+x29==24000];
F=[F,12*x23-x29+x30==24000];
F=[F,12*x24-x30+x31==24000];
F=[F,-8e-5*sqr(x7)+x13==0];
F=[F,-8e-5*sqr(x8)+x14==0];
F=[F,-8e-5*sqr(x9)+x15==0];
F=[F,-8e-5*sqr(x10)+x16==0];
F=[F,-8e-5*sqr(x11)+x17==0];
F=[F,-8e-5*sqr(x12)+x18==0];
F=[F,-4.97*x7+x19==330];
F=[F,-4.97*x8+x20==330];
F=[F,-4.97*x9+x21==330];
F=[F,-4.97*x10+x22==330];
F=[F,-4.97*x11+x23==330];
F=[F,-4.97*x12+x24==330];
F=[F,150<=x1<=1500];
F=[F,150<=x2<=1500];
F=[F,150<=x3<=1500];
F=[F,150<=x4<=1500];
F=[F,150<=x5<=1500];
F=[F,150<=x6<=1500];
F=[F,0<=x7<=1000];
F=[F,0<=x8<=1000];
F=[F,0<=x9<=1000];
F=[F,0<=x10<=1000];
F=[F,0<=x11<=1000];
F=[F,0<=x12<=1000];
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
F=[F,100000==x25];
F=[F,60000<=x26<=120000];
F=[F,60000<=x27<=120000];
F=[F,60000<=x28<=120000];
F=[F,60000<=x29<=120000];
F=[F,60000<=x30<=120000];
F=[F,60000<=x31<=120000];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)-4.3669e+006) <= 1e4)
