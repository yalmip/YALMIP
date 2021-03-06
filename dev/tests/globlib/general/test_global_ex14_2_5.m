function tests = test_global_ex14_2_5
tests = functiontests(localfunctions);

function test1(testCase)
% Model generated from ex14_2_5.gms
% Created 02-Aug-2007 10:38:53 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x5 = sdpvar(1);

% Define objective function 
objective = -(0-x5-0);

% Define constraints 
F = ([]);
F=[F,0.361872516756437*x2/(x1+0.888649896608059*x2)+0.868134622480909*x2/(0.696880695582072*x1+x2)-(0.361872516756437*x1*x2/sqr(x1+0.888649896608059*x2)+0.604986259573375*x2*x1/sqr(0.696880695582072*x1+x2))-2755.64173589155/(219.161+x3)-x5<=-9.20816767045657];
F=[F,0.868134622480909*x1/(0.696880695582072*x1+x2)+0.361872516756437*x1/(x1+0.888649896608059*x2)-(0.321577974600906*x1*x2/sqr(x1+0.888649896608059*x2)+0.868134622480909*x2*x1/sqr(0.696880695582072*x1+x2))-4117.06819797521/(227.438+x3)-x5<=-12.6599269316621];
F=[F,(-0.361872516756437*x2/(x1+0.888649896608059*x2))-0.868134622480909*x2/(0.696880695582072*x1+x2)+0.361872516756437*x1*x2/sqr(x1+0.888649896608059*x2)+0.604986259573375*x2*x1/sqr(0.696880695582072*x1+x2)+2755.64173589155/(219.161+x3)-x5<=9.20816767045657];
F=[F,(-0.868134622480909*x1/(0.696880695582072*x1+x2))-0.361872516756437*x1/(x1+0.888649896608059*x2)+0.321577974600906*x1*x2/sqr(x1+0.888649896608059*x2)+0.868134622480909*x2*x1/sqr(0.696880695582072*x1+x2)+4117.06819797521/(227.438+x3)-x5<=12.6599269316621];
F=[F,x1+x2==1];
F=[F,1e-006<=x1<=1];
F=[F,1e-006<=x2<=1];
F=[F,60<=x3<=100];
F=[F,0<=x5];

% Solve problem
sol = optimize(F,objective,sdpsettings('bmibnb.uppersolver','fmincon','solver','bmibnb'));
testCase.assertTrue(sol.problem==0)
testCase.assertTrue(abs(value(objective)-0) <= 1e-3)