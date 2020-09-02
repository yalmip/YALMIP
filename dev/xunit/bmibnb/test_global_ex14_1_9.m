function tests = test_global_ex14_1_9
tests = functiontests(localfunctions);

function test1(dummy)

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);

% Define objective function 
objective = -(-x2);

% Define constraints 
F = ([]);
F=[F,4510067.11409396*x1*exp(-7548.11926028431/x1)+0.00335570469798658*x1-2020510067.11409*exp(-7548.11926028431/x1)-x2<=1];
F=[F,(-4510067.11409396*x1*exp(-7548.11926028431/x1))-0.00335570469798658*x1+2020510067.11409*exp(-7548.11926028431/x1)-x2<=-1];
F=[F,100<=x1<=1000];

% Gurobi has numericall issues leading to infeasibility during root node
% tightening
sol = optimize(F,objective,sdpsettings('solver','bmibnb','bmibnb.uppersolver','fmincon','allownon',1,'bmibnb.roottight',0));
assert(sol.problem==0)
assert(abs(value(objective)-0) <= 1e-3)