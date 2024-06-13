function tests = test_global_ex14_1_5
tests = functiontests(localfunctions);

function test1(testCase)

yalmip('clear')
sdpvar x1 x2 x3 x4 x5 x6 objvar; 

F = ([]);
F = F + ( - x6 + objvar == 0);
F = F + ( 2*x1 + x2 + x3 + x4 + x5 == 6);
F = F + ( x1 + 2*x2 + x3 + x4 + x5 == 6);
F = F + ( x1 + x2 + 2*x3 + x4 + x5 == 6);
F = F + ( x1 + x2 + x3 + 2*x4 + x5 == 6); 
F = F + ( x1*x2*x3*x4*x5 - x6 <= 1);
F = F + ( - x1*x2*x3*x4*x5 - x6 <= -1);
F = F + ( -2 <= [x1 x2 x3 x4 x5 ] <= 2);

sol = optimize(F,objvar,sdpsettings('bmibnb.uppersolver','fmincon','solver','bmibnb'));

testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(norm(value([x1 x2 x3  x4 x5 ])-[1 1 1 1 1]) <= 1e-5)