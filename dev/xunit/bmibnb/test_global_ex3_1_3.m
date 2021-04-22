function tests = test_global_ex3_1_3
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex3_1_3.gms
% Created 28-Jul-2007 19:00:46 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);
x5 = sdpvar(1);
x6 = sdpvar(1);

objvar = -(-(-25*sqr(x1-2)-sqr(x2-2)-sqr(x3-1)-sqr(x4-4)-sqr(x5-1)-sqr(x6-4)));
% Define constraints 
F = ([]);
%F=[F,-(-25*sqr(x1-2)-sqr(x2-2)-sqr(x3-1)-sqr(x4-4)-sqr(x5-1)-sqr(x6-4))+objvar==0];
F=[F,sqr(x3-3)+x4>=4];
F=[F,sqr(x5-3)+x6>=4];
F=[F,x1-3*x2<=2];
F=[F,-x1+x2<=2];
F=[F,x1+x2<=6];
F=[F,x1+x2>=2];
F=[F,0<=x1];
F=[F,0<=x2];
F=[F,1<=x3<=5];
F=[F,0<=x4<=6];
F=[F,1<=x5<=5];
F=[F,0<=x6<=10];

% Solve problem
sol = optimize(F,objvar,sdpsettings('solver','bmibnb','allownonconvex',1));

assert(sol.problem==0)
assert(abs(value(objvar)- -310) <= 1e-1) 