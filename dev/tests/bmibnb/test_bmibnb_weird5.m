function tests = test_bmibnb_weird5
tests = functiontests(localfunctions);

function test1(testCase)

% Issue #187
yalmip('clear')
sdpvar x y
p = sin(1+y)^2+cos(y*x);
sol = optimize([-1 <= [x y] <= 1, p <= 3],p,sdpsettings('solver','bmibnb'));
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(p)-.5403) <= 1e-2) 