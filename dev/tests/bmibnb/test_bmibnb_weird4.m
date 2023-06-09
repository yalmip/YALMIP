function tests = test_bmibnb_weird4
tests = functiontests(localfunctions);

function test1(testCase)

% Issue #188
yalmip('clear')
sdpvar x y 
p = sin(1+y*x)^2+cos(y*x);
sdpvar z w
sol = optimize([-1 <= [x y z w] <= 1, p <= 3],p,sdpsettings('solver','bmibnb'));

testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(p)-.5403) <= 1e-2) 