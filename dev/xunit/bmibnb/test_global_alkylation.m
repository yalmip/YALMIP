function tests = test_global_abel
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from alkylation.gms
% Created 02-Aug-2007 09:48:12 using YALMIP R20070725

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
objective = -0.063*x4*x7+5.04*x1+0.035*x2+10*x3+3.36*x5+0;

% Define constraints 
F = ([]);
F=[F,x1-1.22*x4+x5==0];
F=[F,-98000*x3/(x4*x9+1000*x3)+x6==0];
F=[F,-(x2+x5)/x1+x8==0];
F=[F,(1.12+0.13167*x8-0.00667*x8^2)*x1-0.99*x4>=0];
F=[F,-(1.12+0.13167*x8-0.00667*x8^2)*x1+1.01010101010101*x4>=0];
F=[F,1.098*x8-0.038*x8^2+0.325*x6-0.99*x7>=-57.425];
F=[F,-(1.098*x8-0.038*x8^2)-0.325*x6+1.01010101010101*x7>=57.425];
F=[F,-0.9*x9-0.222*x10>=-35.82];
F=[F,1.11111111111111*x9+0.222*x10>=35.82];
F=[F,3*x7-0.99*x10>=133];
F=[F,-3*x7+1.01010101010101*x10>=-133];
F=[F,1e-006<=x1<=2000];
F=[F,1e-006<=x2<=16000];
F=[F,1e-006<=x3<=120];
F=[F,1e-006<=x4<=5000];
F=[F,1e-006<=x5<=2000];
F=[F,85<=x6<=93];
F=[F,90<=x7<=95];
F=[F,3<=x8<=12];
F=[F,0.01<=x9<=4];
F=[F,145<=x10<=162];

% Solve problem
%sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1))
%assert(sol.problem == 0)
%assert(abs(value(objective)-  -1.768806963716253e+003) <= 1e-1) 