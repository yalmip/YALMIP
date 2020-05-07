function tests = test_global_ex14_1_3
tests = functiontests(localfunctions);

function test1(dummy)

yalmip('clear');
sdpvar x1 x2 x3 objvar

F = ([]);
F = F + ( - x3 + objvar == 0); 
F = F + ( 10000*x1*x2 - x3 <= 1); 
F = F + ( - 10000*x1*x2 - x3 <= -1); 
F = F + ( exp(-x1) + exp(-x2) - x3 <= 1.001);
F = F + ( (-exp(-x1)) - exp(-x2) - x3 <= -1.001); 
F = F + (5.49e-6 <= x1 <= 4.553) + (18.21 >= x2 >= 0.0021961);

sol = optimize(F,objvar,sdpsettings('solver','bmibnb','bmibnb.absgaptol',1e-8,'bmibnb.relgaptol',1e-8))
assert(sol.problem == 0)
assert(value(objvar) <= 1e-3)
