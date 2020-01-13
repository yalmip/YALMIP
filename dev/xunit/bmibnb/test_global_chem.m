function tests = test_global_chem
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from chem.gms
% Created 21-Aug-2007 17:09:27 using YALMIP R20070810
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
x11 = sdpvar(1);

% Define objective function
objective = -(-(x1*(log(x1/x11)-6.05576803624071)+x2*(log(x2/x11)-17.1307680362407)+x3*(log(x3/x11)-34.0207680362407)+x4*(log(x4/x11)-5.88076803624071)+x5*(log(x5/x11)-24.6877680362407)+x6*(log(x6/x11)-14.9527680362407)+x7*(log(x7/x11)-24.0667680362407)+x8*(log(x8/x11)-10.6747680362407)+x9*(log(x9/x11)-26.6287680362407)+x10*(log(x10/x11)-22.1447680362407))+0-(0));

% Define constraints
F = ([]);
F=[F,x1+2*x2+2*x3+x6+x10==2];
F=[F,x4+2*x5+x6+x7==1];
F=[F,x3+x7+x8+2*x9+x10==1];
F=[F,-x1-x2-x3-x4-x5-x6-x7-x8-x9-x10+x11==0];
F=[F,0.001<=x1];
F=[F,0.001<=x2];
F=[F,0.001<=x3];
F=[F,0.001<=x4];
F=[F,0.001<=x5];
F=[F,0.001<=x6];
F=[F,0.001<=x7];
F=[F,0.001<=x8];
F=[F,0.001<=x9];
F=[F,0.001<=x10];
F=[F,0.01<=x11];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)- -4.771E+001) <=  10)
