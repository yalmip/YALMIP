function tests = test_global_st_qpc_m1
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_qpc-m1.gms
% Created 21-Aug-2007 18:36:51 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);
x5 = sdpvar(1);

% Define objective function 
objective = -(-(10*x1-0.34*x1*x1-0.28*x1*x2+10*x2-0.22*x1*x3+10*x3-0.24*x1*x4+10*x4-0.51*x1*x5+10*x5-0.28*x2*x1-0.34*x2*x2-0.23*x2*x3-0.24*x2*x4-0.45*x2*x5-0.22*x3*x1-0.23*x3*x2-0.35*x3*x3-0.22*x3*x4-0.34*x3*x5-0.24*x4*x1-0.24*x4*x2-0.22*x4*x3-0.2*x4*x4-0.38*x4*x5-0.51*x5*x1-0.45*x5*x2-0.34*x5*x3-0.38*x5*x4-0.99*x5*x5)+0-(0));

% Define constraints 
F = ([]);
F=[F,x1+x2+2*x3+x4+x5>=10];
F=[F,2*x1+3*x2+x5>=8];
F=[F,x2+4*x3-x4+2*x5>=12];
F=[F,8*x1-x2-x3+6*x4>=20];
F=[F,-2*x1-x2-3*x3-x4-x5>=-30];
F=[F,0<=x1];
F=[F,0<=x2];
F=[F,0<=x3];
F=[F,0<=x4];
F=[F,0<=x5];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--4.737777777742748e+002) <= 1e-2) 