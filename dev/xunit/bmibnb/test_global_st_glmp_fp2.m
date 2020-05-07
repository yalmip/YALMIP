function tests = test_global_st_glmp_fp2
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_glmp_fp2.gms
% Created 06-Aug-2007 09:41:33 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);

% Define objective function 
objective = -(-x3*x4+0-(0));

% Define constraints 
F = ([]);
F=[F,2*x1+x2<=14];
F=[F,x1+x2<=10];
F=[F,1.44*x1+x2>=4.89];
F=[F,-1.58*x1+x2<=5.65];
F=[F,-1.03*x1+x2<=5.93];
F=[F,x1+2*x2>=6];
F=[F,x1-x2<=3];
F=[F,x1+x2-x3==0];
F=[F,x1-x2-x4==-7];
F=[F,0<=x1<=5];
F=[F,0<=x2<=10];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)-7.3445) <= 1e-2)