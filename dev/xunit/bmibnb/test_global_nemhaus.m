function tests = test_global_nemhaus
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from nemhaus.gms
% Created 24-Jul-2007 13:45:26 using YALMIP R20070523

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);
x5 = sdpvar(1);
x6 = sdpvar(1);

% Define objective function 
objective = -(-(2*x2*x4+4*x2*x5+3*x2*x6+6*x3*x4+2*x3*x5+3*x3*x6+5*x4*x5+3*x4*x6+3*x5*x6));

% Define constraints 
F = ([]);
F=[F,x2==1];
F=[F,x3==1];
F=[F,x4==1];
F=[F,x5==1];
F=[F,x6==1];
F=[F,0<=x2];
F=[F,0<=x3];
F=[F,0<=x4];
F=[F,0<=x5];
F=[F,0<=x6];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)-31) <= 1e-2)