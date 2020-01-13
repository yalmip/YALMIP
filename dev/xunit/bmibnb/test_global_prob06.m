function tests = test_global_prob06
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from prob06.gms
% Created 17-Mar-2008 10:40:41 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
objvar = sdpvar(1);
x2 = sdpvar(1);

% Define constraints 
F = ([]);
F=[F,0.25*objvar-0.0625*sqr(objvar)-0.0625*sqr(x2)+0.5*x2<=1];
F=[F,0.0714285714285714*sqr(objvar)+0.0714285714285714*sqr(x2)-0.428571428571429*objvar-0.428571428571429*x2<=-1];
F=[F,1<=objvar<=5.5];
F=[F,1<=x2<=5.5];

% Solve problem
sol = optimize(F,objvar,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objvar)-1.1771) <= 1e-1) 