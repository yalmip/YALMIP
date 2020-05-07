function tests = test_global_st_jcbpaf2
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_jcbpaf2.gms
% Created 06-Aug-2007 09:37:26 using YALMIP R20070725

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

% Define objective function 
objective = -(-(x1*x6-x1-x6+x2*x7-2*x2-2*x7+x3*x8-3*x3-3*x8+x4*x9-4*x4-4*x9+x5*x10-5*x5-5*x10)+0-(0));

% Define constraints 
F = ([]);
F=[F,x1+7*x2+5*x3+5*x4-6*x6-3*x7-3*x8+5*x9-7*x10<=80];
F=[F,-3*x1+3*x2+8*x3+7*x4-9*x5-7*x6-9*x7+8*x9-7*x10<=57];
F=[F,x1+x3+3*x4+8*x5+9*x6+9*x8-7*x9-8*x10<=92];
F=[F,-x1-2*x2+2*x3+9*x5+5*x6-3*x7+x8-x9-5*x10<=55];
F=[F,-5*x1+8*x2-8*x3+3*x5+4*x7-5*x8-2*x9+9*x10<=76];
F=[F,4*x1-x2+6*x3-4*x4-7*x5-8*x6-7*x7+6*x8-2*x9-9*x10<=14];
F=[F,7*x2+4*x3+9*x5-6*x8-5*x9-5*x10<=47];
F=[F,-5*x1-x2+7*x4-x5+2*x6+5*x7-8*x8-5*x9+2*x10<=51];
F=[F,-4*x1-7*x2-9*x4+2*x5+6*x6-9*x7+x8-5*x9<=36];
F=[F,-2*x1+6*x2+8*x4-6*x5+8*x6+8*x7+5*x8+2*x9-7*x10<=92];
F=[F,x1+x2+x3-2*x4+x5+x6+x7+4*x8+x9+3*x10<=200];
F=[F,x1+x2+x3+x4+x5>=1];
F=[F,x6+x7+x8+x9+x10>=2];
F=[F,0<=x1<=100];
F=[F,0<=x2<=100];
F=[F,0<=x3<=100];
F=[F,0<=x4<=100];
F=[F,0<=x5<=100];
F=[F,0<=x6<=100];
F=[F,0<=x7<=100];
F=[F,0<=x8<=100];
F=[F,0<=x9<=100];
F=[F,0<=x10<=100];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--794.8559) <= 1e-2)