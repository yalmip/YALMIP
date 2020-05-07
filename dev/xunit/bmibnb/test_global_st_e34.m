function tests = test_global_st_e34
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_e34.gms
% Created 06-Aug-2007 09:45:55 using YALMIP R20070725

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
objective = -(-4.3*x1-31.8*x2-63.3*x3-15.8*x4-68.5*x5-4.7*x6+0-(0));

% Define constraints 
F = ([]);
F=[F,17.1*x1-169*x1*x3+204.2*x3-3580*x3*x5+623.4*x5-3810*x4*x6+212.3*x4+1495.5*x6-18500*x4*x6+38.2*x2>=4.97];
F=[F,17.9*x1-139*x1*x3+113.9*x3-2450*x4*x5+169.7*x4+337.8*x5-16600*x4*x6+1385.2*x6-17200*x5*x6+36.8*x2>=-1.88];
F=[F,26000*x4*x5-70*x4-819*x5-273*x2>=-69.08];
F=[F,159.9*x1-14000*x1*x6+2198*x6-311*x2+587*x4+391*x5>=-118.02];
F=[F,0<=x1<=0.31];
F=[F,0<=x2<=0.046];
F=[F,0<=x3<=0.068];
F=[F,0<=x4<=0.042];
F=[F,0<=x5<=0.028];
F=[F,0<=x6<=0.0134];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)-0.01562) <= 1e-2) 