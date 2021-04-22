function tests = test_global_haverly
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from haverly.gms
% Created 17-Mar-2008 16:35:51 using YALMIP R20070810

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

% Define objective function 
objective = x1-x2-0-(0);

% Define constraints 
F = ([]);
F=[F,x1-6*x3-16*x4-10*x5==0];
F=[F,x2-9*x6-15*x7==0];
F=[F,x6-x8-x10==0];
F=[F,x7-x9-x11==0];
F=[F,x3+x4-x10-x11==0];
F=[F,x5-x8-x9==0];
F=[F,x12*(x10+x11)-3*x3-x4==0];
F=[F,x12*x10-2.5*x10-0.5*x8<=0];
F=[F,x12*x11-1.5*x11+0.5*x9<=0];
F=[F,0<=x1];
F=[F,0<=x2];
F=[F,0<=x3];
F=[F,0<=x4];
F=[F,0<=x5];
F=[F,0<=x6<=100];
F=[F,0<=x7<=200];
F=[F,0<=x8];
F=[F,0<=x9];
F=[F,0<=x10];
F=[F,0<=x11];
F=[F,0<=x12];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1))
assert(sol.problem==0)
assert(abs(value(objective)--400) <= 1e-2) 