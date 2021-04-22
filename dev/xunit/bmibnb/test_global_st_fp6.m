function tests = test_global_st_fp6
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_fp6.gms
% Created 24-Jul-2007 14:16:08 using YALMIP R20070523

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
objective = -(-(48*x1-50*sqr(x1)-50*sqr(x2)+42*x2-50*sqr(x3)+48*x3-50*sqr(x4)+45*x4-50*sqr(x5)+44*x5-50*sqr(x6)+41*x6-50*sqr(x7)+47*x7-50*sqr(x8)+42*x8-50*sqr(x9)+45*x9-50*sqr(x10)+46*x10));

% Define constraints 
F = ([]);
F=[F,-2*x1-6*x2-x3-3*x5-3*x6-2*x7-6*x8-2*x9-2*x10<=-4];
F=[F,6*x1-5*x2+8*x3-3*x4+x6+3*x7+8*x8+9*x9-3*x10<=22];
F=[F,-5*x1+6*x2+5*x3+3*x4+8*x5-8*x6+9*x7+2*x8-9*x10<=-6];
F=[F,9*x1+5*x2-9*x4+x5-8*x6+3*x7-9*x8-9*x9-3*x10<=-23];
F=[F,-8*x1+7*x2-4*x3-5*x4-9*x5+x6-7*x7-x8+3*x9-2*x10<=-12];
F=[F,0<=x1<=1];
F=[F,0<=x2<=1];
F=[F,0<=x3<=1];
F=[F,0<=x4<=1];
F=[F,0<=x5<=1];
F=[F,0<=x6<=1];
F=[F,0<=x7<=1];
F=[F,0<=x8<=1];
F=[F,0<=x9<=1];
F=[F,0<=x10<=1];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--39) <= 1e-2) 