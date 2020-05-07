function tests = test_global_meanvar
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from meanvar.gms
% Created 24-Jul-2007 13:41:24 using YALMIP R20070523

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);
x5 = sdpvar(1);
x6 = sdpvar(1);
x7 = sdpvar(1);
x8 = sdpvar(1);
x9 = sdpvar(1);

% Define objective function 
objective = -(-0.5*(42.18*x3*x3+20.18*x3*x4+10.88*x3*x5+5.3*x3*x6+12.32*x3*x7+23.84*x3*x8+17.41*x3*x9+20.18*x4*x3+70.89*x4*x4+21.58*x4*x5+15.41*x4*x6+23.24*x4*x7+23.8*x4*x8+12.62*x4*x9+10.88*x5*x3+21.58*x5*x4+25.51*x5*x5+9.6*x5*x6+22.63*x5*x7+13.22*x5*x8+4.7*x5*x9+5.3*x6*x3+15.41*x6*x4+9.6*x6*x5+22.33*x6*x6+10.32*x6*x7+10.46*x6*x8+x6*x9+12.32*x7*x3+23.24*x7*x4+22.63*x7*x5+10.32*x7*x6+30.01*x7*x7+16.36*x7*x8+7.2*x7*x9+23.84*x8*x3+23.8*x8*x4+13.22*x8*x5+10.46*x8*x6+16.36*x8*x7+42.23*x8*x8+9.9*x8*x9+17.41*x9*x3+12.62*x9*x4+4.7*x9*x5+x9*x6+7.2*x9*x7+9.9*x9*x8+16.42*x9*x9));

% Define constraints 
F = ([]);
F=[F,x2-0.1287*x3-0.1096*x4-0.0501*x5-0.1524*x6-0.0763*x7-0.1854*x8-0.062*x9==0];
F=[F,x3+x4+x5+x6+x7+x8+x9==1];
F=[F,0.115==x2];
F=[F,0<=x3<=1];
F=[F,0<=x4<=1];
F=[F,0<=x5<=1];
F=[F,0<=x6<=1];
F=[F,0<=x7<=1];
F=[F,0<=x8<=1];
F=[F,0<=x9<=1];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1,'quadprog.Algorithm','interior-point-convex'));
assert(sol.problem==0)
assert(abs(value(objective)-5.2434) <= 1e-2)