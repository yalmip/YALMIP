function tests = test_global_ex7_3_1
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex7_3_1.gms
% Created 17-Mar-2008 11:07:52 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);

% Define objective function 
objective = -(-x4+0-(0));

% Define constraints 
F = ([]);
F=[F,10*sqr(x2)*power(x3,3)+10*power(x2,3)*sqr(x3)+200*sqr(x2)*sqr(x3)+100*power(x2,3)*x3+100*x2*power(x3,3)+x1*x2*sqr(x3)+x1*sqr(x2)*x3+1000*x2*sqr(x3)+8*x1*sqr(x3)+1000*sqr(x2)*x3+8*x1*sqr(x2)+6*x1*x2*x3-sqr(x1)+60*x1*x3+60*x1*x2-200*x1<=0];
F=[F,-x1-800*x4<=-800];
F=[F,x1-800*x4<=800];
F=[F,-x2-2*x4<=-4];
F=[F,x2-2*x4<=4];
F=[F,-x3-3*x4<=-6];
F=[F,x3-3*x4<=6];
F=[F,0<=x1];
F=[F,0<=x2];
F=[F,0<=x3];
F=[F,0<=x4];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1,'bmibnb.uppersolver','fmincon'));
assert(sol.problem==0)
assert(abs(value(objective)-0.3417) <= 1e-2)