function tests = test_global_ex2_1_3
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex2_1_3.gms
% Created 28-Jul-2007 18:43:34 using YALMIP R20070725

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
x11 = sdpvar(1);
x12 = sdpvar(1);
x13 = sdpvar(1);
objvar = sdpvar(1);

% Define constraints 
F = ([]);
F=[F,-(5*x1-0.5*(10*x1*x1+10*x2*x2+10*x3*x3+10*x4*x4)+5*x2+5*x3+5*x4)+x5+x6+x7+x8+x9+x10+x11+x12+x13+objvar==0];
F=[F,2*x1+2*x2+x10+x11<=10];
F=[F,2*x1+2*x3+x10+x12<=10];
F=[F,2*x2+2*x3+x11+x12<=10];
F=[F,-8*x1+x10<=0];
F=[F,-8*x2+x11<=0];
F=[F,-8*x3+x12<=0];
F=[F,-2*x4-x5+x10<=0];
F=[F,-2*x6-x7+x11<=0];
F=[F,-2*x8-x9+x12<=0];
F=[F,0<=x1<=1];
F=[F,0<=x2<=1];
F=[F,0<=x3<=1];
F=[F,0<=x4<=1];
F=[F,0<=x5<=1];
F=[F,0<=x6<=1];
F=[F,0<=x7<=1];
F=[F,0<=x8<=1];
F=[F,0<=x9<=1];
F=[F,0<=x10];
F=[F,0<=x11];
F=[F,0<=x12];
F=[F,0<=x13<=1];

% Solve problem
sol = optimize(F,objvar,sdpsettings('solver','bmibnb','allownonconvex',1));

assert(sol.problem==0)
assert(abs(value(objvar)--15) <= 1e-2)