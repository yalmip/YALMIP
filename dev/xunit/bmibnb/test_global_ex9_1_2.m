function tests = test_global_ex9_1_2
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex9_1_2.gms
% Created 11-Mar-2008 15:58:16 using YALMIP R20070810

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

% Define objective function 
objective = -0-x2-3*x3-(0);

% Define constraints 
F = ([]);
F=[F,-x2+x3+x4==3];
F=[F,x2+2*x3+x5==12];
F=[F,4*x2-x3+x6==12];
F=[F,-x3+x7==0];
F=[F,x8*x4==0];
F=[F,x9*x5==0];
F=[F,x10*x6==0];
F=[F,x11*x7==0];
F=[F,x8+2*x9-x10-x11==-1];
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

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--16) <= 1e-2)