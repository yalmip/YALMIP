function tests = test_bmibnb_weird2
tests = functiontests(localfunctions);

function test1(testCase)

sdpvar x
obj = blackbox(x,@(x)(sin(10*x)+abs(sin(x))+x));
sol = optimize((-pi <= x <= pi),obj,sdpsettings('solver','bmibnb'));

testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(obj)--3.21) <= 2e-1) 

