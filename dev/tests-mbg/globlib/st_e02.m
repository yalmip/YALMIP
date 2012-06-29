% Model generated from st_e02.gms
% Created 21-Aug-2007 18:41:24 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);

% Define objective function 
objective = -(-x3+0-(0));

% Define constraints 
F = set([]);
F=[F,30*x1-6*x1*x1-x3==-250];
F=[F,20*x2-12*x2*x2-x3==-300];
F=[F,0.5*sqr(x1+x2)-x3==-150];
F=[F,0<=x1<=9.422];
F=[F,0<=x2<=5.9023];
F=[F,0<=x3<=267.4171];

% Solve problem
sol = solvesdp(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
mbg_assertfalse(sol.problem)
mbg_asserttolequal(double(objective),2.011593340563166e+002 , 1e-2);