function tests = test_global_alkyl
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from alkyl.gms
% Created 17-Mar-2008 11:04:26 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
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
x15 = sdpvar(1);

% Define objective function 
objective = -(6.3*x5*x8+0-5.04*x2-0.35*x3-x4-3.36*x6-(0));

% Define constraints 
F = ([]);
F=[F,-0.819672131147541*x2+x5-0.819672131147541*x6==0];
F=[F,0.98*x4-x7*(0.01*x5*x10+x4)==0];
F=[F,-x2*x9+10*x3+x6==0];
F=[F,x5*x12-x2*(1.12+0.13167*x9-0.0067*x9*x9)==0];
F=[F,x8*x13-0.01*(1.098*x9-0.038*x9*x9)-0.325*x7==0.57425];
F=[F,x10*x14+22.2*x11==35.82];
F=[F,x11*x15-3*x8==-1.33];
F=[F,0<=x2<=2];
F=[F,0<=x3<=1.6];
F=[F,0<=x4<=1.2];
F=[F,0<=x5<=5];
F=[F,0<=x6<=2];
F=[F,0.85<=x7<=0.93];
F=[F,0.9<=x8<=0.95];
F=[F,3<=x9<=12];
F=[F,1.2<=x10<=4];
F=[F,1.45<=x11<=1.62];
F=[F,0.99<=x12<=1.0101];
F=[F,0.99<=x13<=1.0101];
F=[F,0.9<=x14<=1.1111];
F=[F,0.99<=x15<=1.0101];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--1.7650) <= 1e-2)