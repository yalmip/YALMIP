function tests = test_global_st_rv1
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_rv1.gms
% Created 06-Aug-2007 09:18:22 using YALMIP R20070725

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
x8 = sdpvar(1);
x9 = sdpvar(1);
x10 = sdpvar(1);

% Define objective function 
objective = -(-(-0.00055*sqr(x1)-0.0583*x1-0.0019*sqr(x2)-0.2318*x2-0.0002*sqr(x3)-0.0108*x3-0.00095*sqr(x4)-0.1634*x4-0.0046*sqr(x5)-0.138*x5-0.0035*sqr(x6)-0.357*x6-0.00315*sqr(x7)-0.1953*x7-0.00475*sqr(x8)-0.361*x8-0.0048*sqr(x9)-0.1824*x9-0.003*sqr(x10)-0.162*x10));

% Define constraints 
F = ([]);
F=[F,8*x1+7*x2+9*x3+9*x5+8*x6+2*x7+4*x9+x10<=530];
F=[F,3*x1+4*x2+6*x3+9*x4+6*x6+9*x7+x8+x10<=395];
F=[F,2*x2+x3+5*x4+5*x5+7*x7+4*x8+2*x9<=350];
F=[F,5*x1+7*x3+x4+7*x5+5*x6+7*x8+9*x9+5*x10<=405];
F=[F,x1+x2+x3+x4+x5+x6+x7+x8+x9+x10<=200];
F=[F,0<=x1];
F=[F,0<=x2];
F=[F,0<=x3];
F=[F,0<=x4];
F=[F,0<=x5];
F=[F,0<=x6];
F=[F,0<=x7];
F=[F,0<=x8];
F=[F,0<=x9];
F=[F,0<=x10];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--59.9439) <= 1e-1) 