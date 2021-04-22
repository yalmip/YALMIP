function tests = test_global_st_glmp_ss1
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_glmp_ss1.gms
% Created 06-Aug-2007 09:39:51 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x4 = sdpvar(1);
x5 = sdpvar(1);
x6 = sdpvar(1);

% Define objective function 
objective = -(-x5*x6+0-x4-(0));

% Define constraints 
F = ([]);
F=[F,x1-2*x2<=100];
F=[F,-3*x1-4*x2<=-12];
F=[F,-x1-x2<=100];
F=[F,-x1+4*x2<=100];
F=[F,-x1+2*x2<=18];
F=[F,3*x1+4*x2<=100];
F=[F,x1+x2<=13];
F=[F,x1-4*x2<=8];
F=[F,x1-x4==0];
F=[F,x1-x2-x5==-10];
F=[F,x1+x2-x6==6];
F=[F,0<=x1<=13];
F=[F,0<=x2<=13];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--24.5714) <= 3e-2)