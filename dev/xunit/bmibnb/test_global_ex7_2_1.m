function tests = test_global_ex7_2_1
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex7_2_1.gms
% Created 28-Jul-2007 19:06:49 using YALMIP R20070725

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
objvar = sdpvar(1);

% Define constraints 
F = ([]);
F=[F,-(0.035*x1*x6-0.063*x3*x5+1.715*x1+4.0565*x3)-10*x2+objvar==3000];
F=[F,0.0059553571*sqr(x6)+0.88392857*x3/x1-0.1175625*x6<=1];
F=[F,1.1088*x1/x3+0.1303533*x1/x3*x6-0.0066033*x1/x3*sqr(x6)<=1];
F=[F,0.00066173269*sqr(x6)-0.019120592*x6-0.0056595559*x4+0.017239878*x5<=1];
F=[F,56.85075/x5+1.08702*x6/x5+0.32175*x4/x5-0.03762*sqr(x6)/x5<=1];
F=[F,2462.3121*x2/x3/x4-25.125634*x2/x3+0.006198*x7<=1];
F=[F,161.18996/x7+5000*x2/x3/x7-489510*x2/x3/x4/x7<=1];
F=[F,44.333333/x5+0.33*x7/x5<=1];
F=[F,0.022556*x5-0.007595*x7<=1];
F=[F,-0.0005*x1+0.00061*x3<=1];
F=[F,0.819672*x1/x3+0.819672/x3<=1];
F=[F,24500*x2/x3/x4-250*x2/x3<=1];
F=[F,1.2244898e-5*x3/x2*x4+0.010204082*x4<=1];
F=[F,6.25e-5*x1*x6+6.25e-5*x1-7.625E-5*x3<=1];
F=[F,1.22*x3/x1+1/x1-x6<=1];
F=[F,1500<=x1<=2000];
F=[F,1<=x2<=120];
F=[F,3000<=x3<=3500];
F=[F,85<=x4<=93];
F=[F,90<=x5<=95];
F=[F,3<=x6<=12];
F=[F,145<=x7<=162];

% Solve problem
sol = optimize(F,objvar,sdpsettings('solver','bmibnb','allownonconvex',1));

assert(sol.problem==0 | sol.problem == 3);
assert(abs(value(objvar)-1.12100399602152301e+003) <= 200)