function test_operator_quadratic_in_max

sdpvar x y
obj = -x-y;
F = set(max([x^2+y^2 x+y]) <= 3);
sol = solvesdp(set(max([x^2+y^2 x+y]) <= 3),obj)
assertTrue(sol.problem == 0)
assertElementsAlmostEqual(double(obj),-sqrt(3/2)*2,'absolute', 1e-4);

sdpvar x y
obj = -x-y;
sol = solvesdp(set(0 <= min([-x^2-y^2 -x-y]) +3),obj)
assertTrue(sol.problem == 0)
assertElementsAlmostEqual(double(obj),-sqrt(3/2)*2,'absolute', 1e-4);

sdpvar x y
obj = -x-y;
sol = solvesdp(set(max([1 x y x^2]) <= min([-x^2-y^2 -x-y]) +3),obj)
assertTrue(sol.problem == 0)
assertElementsAlmostEqual(double(obj),-2, 'absolute',1e-4);

sdpvar x y
obj = -x-y;
sol = solvesdp(set(abs([1 x y x^2]) <= min([-x^2-y^2 -x-y]) +3),obj,sdpsettings('allownonconvex',1))
assertTrue(sol.problem == 0)

sdpvar x y
obj = -x-y;
sol = solvesdp(set(abs([1 x y x^2]) <= min([-x^2-y^2 -x-y]) +3),obj,sdpsettings('allownonconvex',0))
assertTrue(sol.problem == 14)



