function tests = test_global_st_e19
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_e19.gms
% Created 22-Aug-2007 09:39:00 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);

% Define objective function 
objective = -(-(power(x1,4)-14*sqr(x1)+24*x1-sqr(x2))+0-(0));

% Define constraints 
F = ([]);
F=[F,-x1+x2<=8];
F=[F,(-sqr(x1))-2*x1+x2<=-2];
F=[F,-8<=x1<=10];
F=[F,0<=x2<=10];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--1.187048598118708e+002) <= 1e-2) 