function tests = test_global_st_e33
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_e33.gms
% Created 06-Aug-2007 09:45:57 using YALMIP R20070725

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

% Define objective function 
objective = -(-6*x1-16*x2+9*x5-10*x6-10*x7+15*x8+0-(0));

% Define constraints 
F = ([]);
F=[F,x1+x2-x3-x4==0];
F=[F,-x9*(x3+x4)+0.03*x1+0.01*x2==0];
F=[F,x3-x5+x6==0];
F=[F,x4+x7-x8==0];
F=[F,x9*x3-0.025*x5+0.02*x6<=0];
F=[F,x9*x4+0.02*x7-0.015*x8<=0];
F=[F,0<=x1<=300];
F=[F,0<=x2<=300];
F=[F,0<=x3<=100];
F=[F,0<=x4<=200];
F=[F,0<=x5<=100];
F=[F,0<=x6<=300];
F=[F,0<=x7<=100];
F=[F,0<=x8<=200];
F=[F,0.01<=x9<=0.03];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--400) <= 1e-2) 