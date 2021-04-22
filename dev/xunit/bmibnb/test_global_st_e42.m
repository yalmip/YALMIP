function tests = test_global_st_e42
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_e42.gms
% Created 06-Aug-2007 09:45:28 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);
x5 = sdpvar(1);
x6 = sdpvar(1);
x7 = sdpvar(1);

% Define objective function 
objective = -(-x1-x2+0-(0));

% Define constraints 
F = ([]);
F=[F,51.5712*x3*x5+20.7324*x3-25.7856*x5+14.9251*x3*x7-22.2988*x7-10.1409*x6*x7+42.3401*x6-x1+x2-6.6425*x4==-72.82];
F=[F,x3==1];
F=[F,0<=x1];
F=[F,0<=x2];
F=[F,0<=x3<=1];
F=[F,0<=x4<=1];
F=[F,-1<=x5<=1];
F=[F,-1<=x6<=1];
F=[F,0<=x7<=3];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)-18.7842) <= 1e-2) 