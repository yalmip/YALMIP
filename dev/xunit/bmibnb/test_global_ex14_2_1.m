function tests = test_global_ex14_2_1
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex14_2_1.gms
% Created 24-Jul-2007 09:56:43 using YALMIP R20070523

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);
x6 = sdpvar(1);

% Define objective function 
objective = -(0-x6-0);

% Define constraints 
F = ([]);
F=[F,log(x1+0.48*x2+0.768*x3)+x1/(x1+0.48*x2+0.768*x3)+1.55*x2/(1.55*x1+x2+0.544*x3)+0.566*x3/(0.566*x1+0.65*x2+x3)+2787.49800065313/(229.664+x4)-x6<=10.7545020354713];
F=[F,log(1.55*x1+x2+0.544*x3)+0.48*x1/(x1+0.48*x2+0.768*x3)+x2/(1.55*x1+x2+0.544*x3)+0.65*x3/(0.566*x1+0.65*x2+x3)+2665.5415812027/(219.726+x4)-x6<=10.6349978691449];
F=[F,log(0.566*x1+0.65*x2+x3)+0.768*x1/(x1+0.48*x2+0.768*x3)+0.544*x2/(1.55*x1+x2+0.544*x3)+x3/(0.566*x1+0.65*x2+x3)+3643.31361767678/(239.726+x4)-x6<=12.9738026256517];
F=[F,(-log(x1+0.48*x2+0.768*x3))-(x1/(x1+0.48*x2+0.768*x3)+1.55*x2/(1.55*x1+x2+0.544*x3)+0.566*x3/(0.566*x1+0.65*x2+x3))-2787.49800065313/(229.664+x4)-x6<=-10.7545020354713];
F=[F,(-log(1.55*x1+x2+0.544*x3))-(0.48*x1/(x1+0.48*x2+0.768*x3)+x2/(1.55*x1+x2+0.544*x3)+0.65*x3/(0.566*x1+0.65*x2+x3))-2665.5415812027/(219.726+x4)-x6<=-10.6349978691449];
F=[F,(-log(0.566*x1+0.65*x2+x3))-(0.768*x1/(x1+0.48*x2+0.768*x3)+0.544*x2/(1.55*x1+x2+0.544*x3)+x3/(0.566*x1+0.65*x2+x3))-3643.31361767678/(239.726+x4)-x6<=-12.9738026256517];
F=[F,x1+x2+x3==1];
F=[F,1e-006<=x1<=1];
F=[F,1e-006<=x2<=1];
F=[F,1e-006<=x3<=1];
F=[F,20<=x4<=80];
F=[F,0<=x6];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)-0) <= 1e-3)