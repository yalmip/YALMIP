function tests = test_global_st_e30
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_e30.gms
% Created 06-Aug-2007 09:45:58 using YALMIP R20070725

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
x11 = sdpvar(1);
x12 = sdpvar(1);
x13 = sdpvar(1);
x14 = sdpvar(1);

% Define objective function 
objective = -(x14+0-(0));

% Define constraints 
F = ([]);
F=[F,-x12*x7-x1+x3==0];
F=[F,-x12*x8-x2+x4==0];
F=[F,(-x13*x7)-x11*x9-x1+x5==0];
F=[F,(-x13*x8)-x11*x10-x2+x6==0];
F=[F,sqr(x7)+sqr(x8)==1];
F=[F,x8+x9==0];
F=[F,-x7+x10==0];
F=[F,-x12+x14<=0];
F=[F,-x11+x14<=0];
F=[F,2*x1+x2>=-1];
F=[F,2*x3+x4>=-1];
F=[F,2*x5+x6>=-1];
F=[F,x1+x2<=1];
F=[F,x3+x4<=1];
F=[F,x5+x6<=1];
F=[F,-1<=x1<=1];
F=[F,-1<=x2<=1];
F=[F,-1<=x3<=1];
F=[F,-1<=x4<=1];
F=[F,-1<=x5<=1];
F=[F,-1<=x6<=1];
F=[F,-1<=x7<=1];
F=[F,-1<=x8<=1];
F=[F,-1<=x9<=1];
F=[F,-1<=x10<=1];
F=[F,0<=x11<=3];
F=[F,0<=x12<=3];
F=[F,0<=x13<=3];
F=[F,0<=x14<=3];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--1.5811) <= 1e-2) 