function test_global_st_jcbpafex
% Model generated from st_jcbpafex.gms
% Created 06-Aug-2007 09:36:40 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);

% Define objective function 
objective = -(-(x1*x2-x1-x2)+0-(0));

% Define constraints 
F = ([]);
F=[F,-6*x1+8*x2<=3];
F=[F,3*x1-x2<=3];
F=[F,0<=x1<=5];
F=[F,0<=x2<=5];

% Solve problem
sol = solvesdp(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
mbg_asserttrue(sol.problem==0)
mbg_asserttolequal(double(objective),-1.0833 , 1e-2);