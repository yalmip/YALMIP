function tests = test_global_st_glmp_kk90
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_glmp_kk90.gms
% Created 06-Aug-2007 09:40:59 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);
x5 = sdpvar(1);

% Define objective function 
objective = -(-x4*x5-x3+0-(0));

% Define constraints 
F = ([]);
F=[F,2*x1+3*x2>=9];
F=[F,3*x1-x2<=8];
F=[F,-x1+2*x2<=8];
F=[F,x1+2*x2<=12];
F=[F,x1-x3==0];
F=[F,x1-x2-x4==-5];
F=[F,x1+x2-x5==1];
F=[F,0<=x1<=12];
F=[F,3<=x2<=6];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)-3) <= 1e-2)