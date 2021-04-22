function tests = test_global_st_e11
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_e11.gms
% Created 21-Aug-2007 18:42:12 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);

% Define objective function 
objective = 35*x1^0.6+35*x2^0.6;

% Define constraints 
F = ([]);
F=[F,600*x1-x1*x3-50*x3==-5000];
F=[F,600*x2+50*x3==15000];
F=[F,0<=x1<=34];
F=[F,0<=x2<=17];
F=[F,0<=x3<=300];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)-189.3116) <= 1e-2) 