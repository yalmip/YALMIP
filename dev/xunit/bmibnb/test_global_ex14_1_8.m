function tests = test_global_ex14_1_8
tests = functiontests(localfunctions);

function test1(dummy)

yalmip('clear')

sdpvar x1 x2 x3 objvar; 

F = ([]);

F = F + ( - x3 + objvar == 0);

F = F + (  (0.0476666666666666 - 0.0649999999999999*x1)*exp(10*x1/(1 + 0.01*x1)) - x1      - x3 <= 0);

F = F + (  x1 - (0.0476666666666666 - 0.0649999999999999*x1)*exp(10*x1/(1 + 0.01*x1))      - x3 <= 0);

F = F + (  (0.143 + (-0.13*x1) - 0.195*x2)*exp(10*x2/(1 + 0.01*x2)) + x1 - 3*x2 - x3  <= 0);

F = F + (  (-(0.143 + (-0.13*x1) - 0.195*x2)*exp(10*x2/(1 + 0.01*x2))) - x1 + 3*x2 - x3 <= 0);

F = F + (-1 <= [x1 x2 ] <= 1);

sol = optimize(F,objvar,sdpsettings('solver','bmibnb','bmibnb.upper','fmincon','allownon',1));
assert(sol.problem==0)
assert(abs(value(objvar)-0) <= 2e-2)