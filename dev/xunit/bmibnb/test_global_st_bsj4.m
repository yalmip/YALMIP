function tests = test_global_st_bsj4
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_bsj4.gms
% Created 17-Mar-2008 11:00:29 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);
x5 = sdpvar(1);
x6 = sdpvar(1);

% Define objective function 
objective = -(-(10.5*x1-1.5*sqr(x1)-sqr(x2)-3.95*x2-sqr(x3)+3*x3-2*sqr(x4)+5*x4-sqr(x5)+1.5*x5-2.5*sqr(x6)-1.5*x6)+0-(0));

% Define constraints 
F = ([]);
F=[F,x1+x2+x3+x4+x5+x6<=500];
F=[F,x1+3*x2+6*x3+2*x4>=50];
F=[F,3*x5+4*x6>=50];
F=[F,x3+2*x4+3*x5+x6<=350];
F=[F,0<=x1<=99];
F=[F,0<=x2<=99];
F=[F,0<=x3<=99];
F=[F,0<=x4<=99];
F=[F,0<=x5<=99];
F=[F,0<=x6<=99];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--70262.05) <= 1e-2) 