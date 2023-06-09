function tests = test_global_bmibnb_sin1
tests = functiontests(localfunctions);

function test1(testCase)

sdpvar x
obj = sin(10*x)+abs(x)
sol = optimize((-pi <= x <= pi),obj,sdpsettings('solver','bmibnb'));

testCase.assertTrue(sol.problem==0)
testCase.assertTrue(abs(value(obj)--0.84792) <= 1e-4)