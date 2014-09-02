function test_global_st_e39
% Model generated from st_e39.gms
% Created 06-Aug-2007 09:45:51 using YALMIP R20070725

% % Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);

% Define objective function 
objective = -(-(-1/(0.1+power(x1-4,2)+power(x2-4,2))-1/(0.2+power(x1-1,2)+power(x2-1,2))-1/(0.2+power(x1-8,2)+power(x2-8,2)))+0-(0));

% Define constraints 
F = ([]);
% Solve problem
sol = solvesdp(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
mbg_assertfalse(sol.problem == 3)