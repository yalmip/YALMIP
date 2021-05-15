function tests = test_operator_quadrtic_in_max
tests = functiontests(localfunctions);

function test1(testCase)
sdpvar x y
obj = -x-y;
F = (max([x^2+y^2 x+y]) <= 3);
sol = optimize((max([x^2+y^2 x+y]) <= 3),obj)
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(obj)--sqrt(3/2)*2) <= 1e-4);

obj = -x-y;
sol = optimize((0 <= min([-x^2-y^2 -x-y]) +3),obj)
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(obj)--sqrt(3/2)*2) <= 1e-4);

obj = -x-y;
sol = optimize((max([1 x y x^2]) <= min([-x^2-y^2 -x-y]) +3),obj)
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(obj)--2) <= 1e-4);

obj = -x-y;
sol = optimize((abs([1 x y x^2]) <= min([-x^2-y^2 -x-y]) +3),obj,sdpsettings('allownonconvex',1))
testCase.assertTrue(sol.problem == 0)

obj = -x-y;
sol = optimize((abs([1 x y x^2]) <= min([-x^2-y^2 -x-y]) +3),obj,sdpsettings('allownonconvex',0))
testCase.assertTrue(sol.problem == 14)