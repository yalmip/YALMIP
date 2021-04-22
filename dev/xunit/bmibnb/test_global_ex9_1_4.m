function tests = test_global_ex9_1_4
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex9_1_4.gms
% Created 11-Mar-2008 15:55:48 using YALMIP R20070810

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
objective = -0+x2-4*x3-(0);

% Define constraints 
F = ([]);
F=[F,-2*x2+x3+x4==0];
F=[F,2*x2+5*x3+x5==108];
F=[F,2*x2-3*x3+x6==-4];
F=[F,-x3+x7==0];
F=[F,x8*x4==0];
F=[F,x9*x5==0];
F=[F,x10*x6==0];
F=[F,x11*x7==0];
F=[F,x8+5*x9-3*x10-x11==-1];
F=[F,0<=x2];
F=[F,0<=x3];
F=[F,0<=x4<=200];
F=[F,0<=x5<=200];
F=[F,0<=x6<=200];
F=[F,0<=x7<=200];
F=[F,0<=x8<=200];
F=[F,0<=x9<=200];
F=[F,0<=x10<=200];
F=[F,0<=x11<=200];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--37) <= 1e-2)