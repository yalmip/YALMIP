function tests = test_global_ex3_1_1
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex3_1_1.gms
% Created 28-Jul-2007 18:58:50 using YALMIP R20070725

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
objvar = sdpvar(1);

objvar = x1+x2+x3;
% Define constraints 
F = ([]);
%F=[F,-x1-x2-x3+objvar==0];
F=[F,0.0025*x4+0.0025*x6<=1];
F=[F,-0.0025*x4+0.0025*x5+0.0025*x7<=1];
F=[F,-0.01*x5+0.01*x8<=1];
F=[F,100*x1-x1*x6+833.33252*x4<=83333.333];
F=[F,x2*x4-x2*x7-1250*x4+1250*x5<=0];
F=[F,x3*x5-x3*x8-2500*x5<=-1250000];
F=[F,100<=x1<=10000];
F=[F,1000<=x2<=10000];
F=[F,1000<=x3<=10000];
F=[F,10<=x4<=1000];
F=[F,10<=x5<=1000];
F=[F,10<=x6<=1000];
F=[F,10<=x7<=1000];
F=[F,10<=x8<=1000];

% Solve problem
sol = optimize(F,objvar,sdpsettings('bmibnb.upper','fmincon','solver','bmibnb','allownonconvex',1));

assert(sol.problem==0)
assert(abs(value(objvar)- 7049.248) <=  .1)