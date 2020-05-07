function tests = test_global_ex14_1_4
tests = functiontests(localfunctions);

function test1(dummy)

yalmip('clear');
sdpvar x1 x2 x3 objvar

F = ([]);

F = F + (0.5*sin(x1*x2) - 0.5*x1 - 0.0795774703703634*x2 - x3 <= 0);
F = F + (0.920422529629637*exp(2*x1) - 5.4365636*x1 + 0.865255957591193*x2 - x3 <= 2.5019678106022);
F = F + (0.5*x1 - 0.5*sin(x1*x2) + 0.0795774703703634*x2 - x3 <= 0);
F = F + (- x3 + objvar == 0);
F = F + (5.4365636*x1 - 0.920422529629637*exp(2*x1) - 0.865255957591193*x2 - x3 <= -2.5019678106022); 

F = F + (x1>= 0.25) + (x1<= 1)+(6.28 >= x2 >= 1.5);

sol = optimize(F,x3,sdpsettings('solver','bmibnb','bmibnb.upper','fmincon','allownon',1,'bmibnb.absgaptol',1e-8,'bmibnb.relgaptol',1e-8));
assert(sol.problem == 0)
assert(abs(value(objvar)) <= 1e-4)
