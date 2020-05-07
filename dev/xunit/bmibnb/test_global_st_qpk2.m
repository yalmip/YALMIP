function tests = test_global_st_qpk2
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_qpk2.gms
% Created 28-Jul-2007 18:16:25 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);
x5 = sdpvar(1);
x6 = sdpvar(1);

% Define objective function 
objective = -(-(0.5*x1*x2-x1*x1+0.5*x2*x1-x2*x2+0.5*x2*x3+0.5*x3*x2-x3*x3+0.5*x3*x4+0.5*x4*x3-x4*x4+0.5*x4*x5+0.5*x5*x4-x5*x5+0.5*x5*x6+0.5*x6*x5-x6*x6));

% Define constraints 
F = ([]);
F=[F,-x1-2*x2-3*x3-4*x4-5*x5-6*x6<=0];
F=[F,-2*x1-3*x2-4*x3-5*x4-6*x5-x6<=0];
F=[F,-3*x1-4*x2-5*x3-6*x4-x5-2*x6<=0];
F=[F,-4*x1-5*x2-6*x3-x4-2*x5-3*x6<=0];
F=[F,-5*x1-6*x2-x3-2*x4-3*x5-4*x6<=0];
F=[F,-6*x1-x2-2*x3-3*x4-4*x5-5*x6<=0];
F=[F,x1+2*x2+3*x3+4*x4+5*x5+6*x6<=21];
F=[F,2*x1+3*x2+4*x3+5*x4+6*x5+x6<=21];
F=[F,3*x1+4*x2+5*x3+6*x4+x5+2*x6<=21];
F=[F,4*x1+5*x2+6*x3+x4+2*x5+3*x6<=21];
F=[F,5*x1+6*x2+x3+2*x4+3*x5+4*x6<=21];
F=[F,6*x1+x2+2*x3+3*x4+4*x5+5*x6<=21];
F=[F,0<=x1];
F=[F,0<=x2];
F=[F,0<=x3];
F=[F,0<=x4];
F=[F,0<=x5];
F=[F,0<=x6];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));

assert(sol.problem==0)
assert(abs(value(objective)--12.25) <= 1e-2) 