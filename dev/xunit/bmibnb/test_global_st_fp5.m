function tests = test_global_st_fp5
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_fp5.gms
% Created 24-Jul-2007 14:16:06 using YALMIP R20070523

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
objective = -(-(-5*sqr(x1)-20*x1-5*sqr(x2)-80*x2-5*sqr(x3)-20*x3-5*sqr(x4)-50*x4-5*sqr(x5)-60*x5-5*sqr(x6)-90*x6-5*sqr(x7))-10*x8-10*x9-10*x10);

% Define constraints 
F = ([]);
F=[F,-2*x1-6*x2-x3-3*x5-3*x6-2*x7-6*x8-2*x9-2*x10<=-4];
F=[F,6*x1-5*x2+8*x3-3*x4+x6+3*x7+8*x8+9*x9-3*x10<=22];
F=[F,-5*x1+6*x2+5*x3+3*x4+8*x5-8*x6+9*x7+2*x8-9*x10<=-6];
F=[F,9*x1+5*x2-9*x4+x5-8*x6+3*x7-9*x8-9*x9-3*x10<=-23];
F=[F,-8*x1+7*x2-4*x3-5*x4-9*x5+x6-7*x7-x8+3*x9-2*x10<=-12];
F=[F,-7*x1-5*x2-2*x3-6*x5-6*x6-7*x7-6*x8+7*x9+7*x10<=-3];
F=[F,x1-3*x2-3*x3-4*x4-x5-4*x7+x8+6*x9<=1];
F=[F,x1-2*x2+6*x3+9*x4-7*x6+9*x7-9*x8-6*x9+4*x10<=12];
F=[F,-4*x1+6*x2+7*x3+2*x4+2*x5+6*x7+6*x8-7*x9+4*x10<=15];
F=[F,x1+x2+x3+x4+x5+x6+x7+x8+x9+x10<=9];
F=[F,-x1-x2-x3-x4-x5-x6-x7-x8-x9-x10<=-1];
F=[F,0<=x1<=1];
F=[F,0<=x2<=1];
F=[F,0<=x3<=1];
F=[F,0<=x4<=1];
F=[F,0<=x5<=1];
F=[F,0<=x6<=1];
F=[F,0<=x7<=1];
F=[F,0<=x8<=1];
F=[F,0<=x9<=1];
F=[F,0<=x10<=1];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--268.0146) <= 1e-2) 