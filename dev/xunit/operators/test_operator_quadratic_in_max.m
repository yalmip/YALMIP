function tests = test_operator_quadrtic_in_max
tests = functiontests(localfunctions);

function test1(dummy)
sdpvar x y
obj = -x-y;
F = (max([x^2+y^2 x+y]) <= 3);
sol = solvesdp((max([x^2+y^2 x+y]) <= 3),obj)
assert(sol.problem == 0)
assert(abs(double(obj)--sqrt(3/2)*2) <= 1e-4);

sdpvar x y
obj = -x-y;
sol = solvesdp((0 <= min([-x^2-y^2 -x-y]) +3),obj)
assert(sol.problem == 0)
assert(abs(double(obj)--sqrt(3/2)*2) <= 1e-4);

sdpvar x y
obj = -x-y;
sol = solvesdp((max([1 x y x^2]) <= min([-x^2-y^2 -x-y]) +3),obj)
assert(sol.problem == 0)
assert(abs(double(obj)--2) <= 1e-4);

sdpvar x y
obj = -x-y;
sol = solvesdp((abs([1 x y x^2]) <= min([-x^2-y^2 -x-y]) +3),obj,sdpsettings('allownonconvex',1))
assert(sol.problem == 0)

sdpvar x y
obj = -x-y;
sol = solvesdp((abs([1 x y x^2]) <= min([-x^2-y^2 -x-y]) +3),obj,sdpsettings('allownonconvex',0))
assert(sol.problem == 14)



