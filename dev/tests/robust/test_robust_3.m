function tests = test_robust_3
tests = functiontests(localfunctions);

function test1(testCase)
sdpvar x u

sol = optimize(((1+u)*x+x>=0.2) + (0.1 <= u <= 0.3) + (uncertain(u)),x,sdpsettings('verbose',0));

testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(x) - 9.523809524e-002) <= 1e-5)
