function tests = test_global_st_e04
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_e04.gms
% Created 21-Aug-2007 18:41:57 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);

% Define objective function 
objective = -(-(400*x1^0.9+22*(x2-14.7)^1.2)-x3+0-(1000));

% Define constraints 
F = ([]);
F=[F,x3*x1+144*x4>=11520];
F=[F,-exp(11.86-3950/(460+x4))+x2==0];
F=[F,0<=x1<=15.1];
F=[F,14.7<=x2<=94.2];
F=[F,0<=x3<=5371];
F=[F,-459.67<=x4<=80];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)-5.194866244203778e+003) <= 1e-1) 