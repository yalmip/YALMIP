function test_global_st_e09
% Model generated from st_e09.gms
% Created 21-Aug-2007 18:42:08 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);

% Define objective function 
objective = -(2*x1*x2+0-(0));

% Define constraints 
F = ([]);
F=[F,4*x1*x2+2*x1+2*x2<=3];
F=[F,0<=x1<=1];
F=[F,0<=x2<=1];

% Solve problem
sol = solvesdp(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
mbg_asserttrue(sol.problem==0)
mbg_asserttolequal(double(objective),-0.5 , 1e-2);