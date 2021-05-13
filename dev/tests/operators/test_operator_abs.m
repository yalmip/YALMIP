function tests = test_operator_abs
tests = functiontests(localfunctions);

function test1(testCase)

sdpvar x y
obj = abs(1+abs(x-5));
sol = optimize([],obj)

testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(obj) - 1) <= 1e-4);

function test2(testCase)
sdpvar x y
F = (abs(1+abs(x-5)) + abs(y)<=10) 
obj = -x
sol = optimize(F,obj)
testCase.assertTrue(abs(value(obj)--14) <= 1e-4);

