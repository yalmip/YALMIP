function tests = test_global_st_e03
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_e03.gms
% Created 21-Aug-2007 18:41:56 using YALMIP R20070810

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
objective = -(0.063*x4*x7-5.04*x1-0.035*x2-10*x3-3.36*x5+0-(0));

% Define constraints 
F = ([]);
F=[F,x1-1.22*x4+x5==0];
F=[F,x9+0.222*x10==35.82];
F=[F,3*x7-x10==133];
F=[F,0.038*sqr(x8)-1.098*x8-0.325*x6+x7==57.425];
F=[F,x4*x9*x6+1000*x3*x6-98000*x3==0];
F=[F,-x1*x8+x2+x5==0];
F=[F,0.13167*x8*x1+1.12*x1-0.00667*sqr(x8)*x1-x4>=0];
F=[F,1<=x1<=2000];
F=[F,1<=x2<=16000];
F=[F,0<=x3<=120];
F=[F,1<=x4<=5000];
F=[F,0<=x5<=2000];
F=[F,85<=x6<=93];
F=[F,90<=x7<=95];
F=[F,3<=x8<=12];
F=[F,1.2<=x9<=4];
F=[F,145<=x10<=162];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--1.161336609382943e+003) <= 1e-2)