function tests = test_global_circle
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from circle.gms
% Created 17-Mar-2008 11:05:04 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
objvar = sdpvar(1);

% Define constraints 
F = ([]);
F=[F,sqr(2.545724188-x1)+sqr(9.983058643-x2)-sqr(objvar)<=0];
F=[F,sqr(8.589400372-x1)+sqr(6.208600402-x2)-sqr(objvar)<=0];
F=[F,sqr(5.953378204-x1)+sqr(9.920197351-x2)-sqr(objvar)<=0];
F=[F,sqr(3.710241136-x1)+sqr(7.860254203-x2)-sqr(objvar)<=0];
F=[F,sqr(3.629909053-x1)+sqr(2.176232347-x2)-sqr(objvar)<=0];
F=[F,sqr(3.016475803-x1)+sqr(6.757468831-x2)-sqr(objvar)<=0];
F=[F,sqr(4.148474536-x1)+sqr(2.435660776-x2)-sqr(objvar)<=0];
F=[F,sqr(8.706433123-x1)+sqr(3.250724797-x2)-sqr(objvar)<=0];
F=[F,sqr(1.604023507-x1)+sqr(7.020357481-x2)-sqr(objvar)<=0];
F=[F,sqr(5.501896021-x1)+sqr(4.918207429-x2)-sqr(objvar)<=0];
F=[F,0<=objvar];

% Solve problem
sol = optimize(F+(recover(depends(F))<=100),objvar,sdpsettings('solver','bmibnb','allownonconvex',1,'bmibnb.absgaptol',1e-3,'bmibnb.relgaptol',1e-3))
assert(sol.problem==0)
assert(abs(value(objvar)-4.5742) <= 1e-2)