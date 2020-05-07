function tests = test_global_ex3_1_2
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex3_1_2.gms
% Created 28-Jul-2007 18:59:55 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);
x5 = sdpvar(1);
objvar = sdpvar(1);

% Define constraints 
objvar = -40792.141+(0.8356891*x1*x5+37.293239*x1+5.3578547*x3*x3);
F = ([]);
F=[F,0.0056858*x2*x5-0.0022053*x3*x5+0.0006262*x1*x4<=6.665593];
F=[F,0.0022053*x3*x5-0.0056858*x2*x5-0.0006262*x1*x4<=85.334407];
F=[F,0.0071317*x2*x5+0.0021813*x3*x3+0.0029955*x1*x2<=29.48751];
F=[F,(-0.0071317*x2*x5)-0.0021813*x3*x3-0.0029955*x1*x2<=-9.48751];
F=[F,0.0047026*x3*x5+0.0019085*x3*x4+0.0012547*x1*x3<=15.599039];
F=[F,(-0.0047026*x3*x5)-0.0019085*x3*x4-0.0012547*x1*x3<=-10.699039];
F=[F,78<=x1<=102];
F=[F,33<=x2<=45];
F=[F,27<=x3<=45];
F=[F,27<=x4<=45];
F=[F,27<=x5<=45];

% Solve problem
sol = optimize(F,objvar,sdpsettings('solver','bmibnb','bmibnb.upper','fmincon','allownonconvex',1));

assert(sol.problem==0)
assert(abs(value(objvar)--30665.5387) <=  10)