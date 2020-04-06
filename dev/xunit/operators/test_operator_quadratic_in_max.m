function tests = test_operator_quadrtic_in_max
tests = functiontests(localfunctions);

function test1(dummy)
sdpvar x y
obj = -x-y;
F = (max([x^2+y^2 x+y]) <= 3);
sol = optimize((max([x^2+y^2 x+y]) <= 3),obj)
assert(sol.problem == 0)
assert(abs(value(obj)--sqrt(3/2)*2) <= 1e-4);

function test2(dummy)
sdpvar x y
obj = -x-y;
sol = optimize((0 <= min([-x^2-y^2 -x-y]) +3),obj)
assert(sol.problem == 0)
assert(abs(value(obj)--sqrt(3/2)*2) <= 1e-4);

function test3(dummy)
sdpvar x y
obj = -x-y;
sol = optimize((max([1 x y x^2]) <= min([-x^2-y^2 -x-y]) +3),obj)
assert(sol.problem == 0)
assert(abs(value(obj)--2) <= 1e-4);

function test4(dummy)
sdpvar x y
obj = -x-y;
sol = optimize((abs([1 x y x^2]) <= min([-x^2-y^2 -x-y]) +3),obj,sdpsettings('allownonconvex',1))
assert(sol.problem == 0)

function test5(dummy)
sdpvar x y
obj = -x-y;
sol = optimize((abs([1 x y x^2]) <= min([-x^2-y^2 -x-y]) +3),obj,sdpsettings('allownonconvex',0))
assert(sol.problem == 14)