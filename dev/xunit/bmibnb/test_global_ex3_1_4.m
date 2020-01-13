function tests = test_global_ex3_1_4
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex3_1_4.gms
% Created 28-Jul-2007 19:01:03 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
objvar = sdpvar(1);

% Define constraints 
F = ([]);
F=[F,2*x1-x2+x3+objvar==0];
F=[F,x1*(4*x1-2*x2+2*x3)+x2*(2*x2-2*x1-x3)+x3*(2*x1-x2+2*x3)-20*x1+9*x2-13*x3>=-24];
F=[F,x1+x2+x3<=4];
F=[F,3*x2+x3<=6];
F=[F,0<=x1<=2];
F=[F,0<=x2];
F=[F,0<=x3<=3];

% Solve problem
sol = optimize(F,objvar,sdpsettings('solver','bmibnb','allownonconvex',1));

assert(sol.problem==0)
assert(abs(value(objvar)- -4) <= 1e-1)