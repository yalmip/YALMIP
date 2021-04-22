function tests = test_global_st_e41
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_e41.gms
% Created 06-Aug-2007 09:45:49 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);

% Define objective function 
objective = -(-(200*x1^0.6+200*x2^0.6+200*x3^0.6+300*x4^0.6)+0-(0));

% Define constraints 
F = ([]);
F=[F,-(-x3*sqr(1-x1)*sqr(1-x4)-(1-x3)*sqr(1-x2*(1-(1-x1)*(1-x4))))<=0.1];
F=[F,(-x3*sqr(1-x1)*sqr(1-x4))-(1-x3)*sqr(1-x2*(1-(1-x1)*(1-x4)))<=0];
F=[F,0.5<=x1<=1];
F=[F,0.5<=x2<=1];
F=[F,0.5<=x3<=1];
F=[F,0.5<=x4<=1];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)-641.8236) <= 10)