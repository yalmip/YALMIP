function ex14_1_9

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);

% Define objective function 
objective = -(-x2);

% Define constraints 
F = ([]);
F=[F,4510067.11409396*x1*exp(-7548.11926028431/x1)+0.00335570469798658*x1-2020510067.11409*exp(-7548.11926028431/x1)-x2<=1];
F=[F,(-4510067.11409396*x1*exp(-7548.11926028431/x1))-0.00335570469798658*x1+2020510067.11409*exp(-7548.11926028431/x1)-x2<=-1];
F=[F,100<=x1<=1000];


sol = solvesdp(F,objective,sdpsettings('solver','bmibnb','bmibnb.upper','fmincon','allownon',1));
mbg_asserttrue(sol.problem==0);
mbg_asserttolequal(double(objective), 0, 1e-3);