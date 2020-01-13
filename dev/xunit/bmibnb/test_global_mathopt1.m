function tests = test_global_mathopt1
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from mathopt1.gms
% Created 24-Jul-2007 13:38:40 using YALMIP R20070523

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);

% Define objective function 
objective = -(-(10*sqr(sqr(x1)-x2)+sqr((-1)+x1)));

% Define constraints 
F = ([]);
F=[F,x1-x1*x2==0];
F=[F,3*x1+4*x2<=25];
F=[F,-10<=x1<=20];
F=[F,-15<=x2<=20];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)-0) <= 1e-2)