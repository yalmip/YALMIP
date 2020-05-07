function tests = test_global_st_e20
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_e20.gms
% Created 22-Aug-2007 09:36:44 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);
x5 = sdpvar(1);
x6 = sdpvar(1);

% Define objective function 
objective = -(x4+0-(0));

% Define constraints 
F = ([]);
F=[F,0.09755988*x1*x5+x1==1];
F=[F,0.0965842812*x2*x6+x2-x1==0];
F=[F,0.0391908*x3*x5+x3+x1==1];
F=[F,0.03527172*x4*x6+x4-x1+x2-x3==0];
F=[F,x5^0.5+x6^0.5<=4];
F=[F,0<=x1<=1];
F=[F,0<=x2<=1];
F=[F,0<=x3<=1];
F=[F,0<=x4<=1];
F=[F,1e-005<=x5<=16];
F=[F,1e-005<=x6<=16];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--0.388781907395925) <= 1e-2)