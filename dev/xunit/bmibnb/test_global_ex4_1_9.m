function tests = test_global_ex4_1_9
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex4_1_9.gms
% Created 28-Jul-2007 18:55:16 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
objvar = sdpvar(1);

% Define constraints 
objvar = -x1-x2;
F = ([]);
%F=[F,x1+x2+objvar==0];
F=[F,8*power(x1,3)-2*power(x1,4)-8*sqr(x1)+x2<=2];
F=[F,32*power(x1,3)-4*power(x1,4)-88*sqr(x1)+96*x1+x2<=36];
F=[F,0<=x1<=3];
F=[F,0<=x2<=4];

% Solve problem
sol = optimize(F,objvar,sdpsettings('solver','bmibnb','allownonconvex',1));

assert(sol.problem==0)
assert(abs(value(objvar)--5.508) <= 1e-2)