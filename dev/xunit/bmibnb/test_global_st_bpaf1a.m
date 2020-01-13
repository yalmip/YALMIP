function tests = test_global_st_bpaf1a
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_bpaf1a.gms
% Created 17-Mar-2008 10:58:01 using YALMIP R20070810

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
objective = -(-(x1*x6+2*x1+3*x6+x2*x7-4*x2-x7+x3*x8+8*x3-2*x8+x4*x9+4*x4-4*x9+x5*x10+9*x5+5*x10)+0-(0));

% Define constraints 
F = ([]);
F=[F,-8*x1-6*x3+7*x4-7*x5<=1];
F=[F,-6*x1+2*x2-3*x3+9*x4-3*x5<=3];
F=[F,6*x1-7*x3-8*x4+2*x5<=5];
F=[F,-x1+x2-8*x3-7*x4-5*x5<=4];
F=[F,4*x1-7*x2+4*x3+5*x4+x5<=0];
F=[F,5*x7-4*x8+9*x9-7*x10<=0];
F=[F,7*x6+4*x7+3*x8+7*x9+5*x10<=7];
F=[F,6*x6+x7-8*x8+8*x9<=3];
F=[F,-3*x6+2*x7+7*x8+x10<=6];
F=[F,-2*x6-3*x7+8*x8+5*x9-2*x10<=2];
F=[F,0<=x1<=20];
F=[F,0<=x2<=20];
F=[F,0<=x3<=20];
F=[F,0<=x4<=20];
F=[F,0<=x5<=20];
F=[F,0<=x6<=20];
F=[F,0<=x7<=20];
F=[F,0<=x8<=20];
F=[F,0<=x9<=20];
F=[F,0<=x10<=20];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--45.3797) <= 1e-2)