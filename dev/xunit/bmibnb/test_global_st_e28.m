function tests = test_global_st_e28
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_e28.gms
% Created 24-Jul-2007 14:13:39 using YALMIP R20070523

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

% Define objective function 
objective = -(-(5.3578547*sqr(x7)+0.8356891*x5*x9+37.293239*x5)-5000*x4+0-(-40792.141));

% Define constraints 
F = ([]);
F=[F,5*x4-x5+7*x7-x9>=0];
F=[F,-(0.0056858*x6*x9+0.0006262*x5*x8-0.0022053*x7*x9)+x1+2*x4==85.334407];
F=[F,-(0.0071317*x6*x9+0.0029955*x5*x6+0.0021813*sqr(x7))+x2==80.51249];
F=[F,-(0.0047026*x7*x9+0.0012547*x5*x7+0.0019085*x7*x8)+x3+4*x4==9.300961];
F=[F,0<=x1<=92];
F=[F,90<=x2<=110];
F=[F,20<=x3<=25];
F=[F,0<=x4];
F=[F,78<=x5<=102];
F=[F,33<=x6<=45];
F=[F,27<=x7<=45];
F=[F,27<=x8<=45];
F=[F,27<=x9<=45];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--3.0666e+004) <= 10)