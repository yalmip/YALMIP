function tests = test_global_ex9_1_8
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex9_1_8.gms
% Created 24-Jul-2007 12:44:24 using YALMIP R20070523

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
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

% Define objective function 
objective = -0-2*x2+x3+0.5*x4-0;

% Define constraints 
F = ([]);
F=[F,x2+x3<=2];
F=[F,-2*x2+x4-x5+x6==-2.5];
F=[F,x2-3*x3+x5+x7==2];
F=[F,-x4+x8==0];
F=[F,-x5+x9==0];
F=[F,x11*x6==0];
F=[F,x12*x7==0];
F=[F,x13*x8==0];
F=[F,x14*x9==0];
F=[F,x15*x10==0];
F=[F,x11-x13==4];
F=[F,x11+x12-x14==-1];
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

% Solve problem
optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
sol = optimize(F+(-100<=recover(depends(F))<=100),objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--3.25) <= 1e-2)