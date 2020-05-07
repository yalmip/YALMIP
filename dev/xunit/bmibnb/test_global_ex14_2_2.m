function tests = test_global_ex14_2_2
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex14_2_2.gms
% Created 02-Aug-2007 10:33:36 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x5 = sdpvar(1);

% Define objective function 
objective = -(0-x5-0);

% Define constraints 
F = ([]);
F=[F,log(x1+0.191987347447993*x2)+x1/(x1+0.191987347447993*x2)+0.315693799947296*x2/(0.315693799947296*x1+x2)+3643.31361767678/(239.726+x3)-x5<=12.9738026256517];
F=[F,log(0.315693799947296*x1+x2)+0.191987347447993*x1/(x1+0.191987347447993*x2)+x2/(0.315693799947296*x1+x2)+2755.64173589155/(219.161+x3)-x5<=10.2081676704566];
F=[F,(-log(x1+0.191987347447993*x2))-(x1/(x1+0.191987347447993*x2)+0.315693799947296*x2/(0.315693799947296*x1+x2))-3643.31361767678/(239.726+x3)-x5<=-12.9738026256517];
F=[F,(-log(0.315693799947296*x1+x2))-(0.191987347447993*x1/(x1+0.191987347447993*x2)+x2/(0.315693799947296*x1+x2))-2755.64173589155/(219.161+x3)-x5<=-10.2081676704566];
F=[F,x1+x2==1];
F=[F,1e-006<=x1<=1];
F=[F,1e-006<=x2<=1];
F=[F,20<=x3<=80];
F=[F,0<=x5];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)-0) <= 1e-2)