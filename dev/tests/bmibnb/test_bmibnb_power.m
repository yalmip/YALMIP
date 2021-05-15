function tests = test_global_power
tests = functiontests(localfunctions);

function test1(testCase)

sdpvar x y
ops =sdpsettings('solver','bmibnb');
sol = optimize([0 <= [x] <= 5, -5 <= y <= 5, x==.1], (x^y-10)^2,ops)

testCase.assertTrue(sol.problem==0)
testCase.assertTrue(abs(value(y)--1) <= 1e-3)

function test2(testCase)

x = sdpvar(2,1);
y = sdpvar(2,1);
ops =sdpsettings('solver','bmibnb');
sol = optimize([0 <= [x] <= 5, -5 <= y <= 5, x==.1], sum(sum((x.^y-10).^2)),ops)

testCase.assertTrue(sol.problem==0)
testCase.assertTrue(norm(value(y)-[-1;-1]) <= 1e-3)