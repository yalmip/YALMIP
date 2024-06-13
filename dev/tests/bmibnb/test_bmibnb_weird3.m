function tests = test_bmibnb_weird3
tests = functiontests(localfunctions);

function test1(testCase)

sdpvar x
obj = sin(sin(x.^2) + x.^2)+0.01*x.^2-sin(x);
sol = optimize([-2*pi <= x <= 2*pi],obj,sdpsettings('solver','bmibnb'))

testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(obj)--1.734) <= 1e-3)