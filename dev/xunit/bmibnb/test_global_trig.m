function test_global_trig
% Model generated from trig.gms
% Created 28-Jul-2007 18:14:21 using YALMIP R20070725

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
objvar = sdpvar(1);

% Define constraints 
F = ([]);
objvar = (sin(11*x1)+cos(13*x1)-sin(17*x1)-cos(19*x1));
F=[F,5*sin(x1)-x1<=0];
F=[F,-2<=x1<=5];

% Solve problem
sol = solvesdp(F,objvar,sdpsettings('solver','bmibnb','allownonconvex',1));

mbg_asserttrue(sol.problem==0);
mbg_asserttolequal(double(objvar),-3.7625, 1e-2);