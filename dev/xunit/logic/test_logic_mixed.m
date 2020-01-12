function tests = test_logic_mixed
tests = functiontests(localfunctions);

function test1(dummy)
binvar d1 d2

sol =  solvesdp(~(d1==d2),d2)
assert(sol.problem == 0);
assert(abs(double(d1 + d2)-1) <= 1e-3)

sol =  solvesdp((d1~=d2),d1+d2)
assert(sol.problem == 0);
assert(abs(double(d1 + d2)-1) <= 1e-3)

sol =  solvesdp(iff(d1,d2),d1+d2)
assert(sol.problem == 0);
assert(abs(double(d1 + d2)-0) <= 1e-3)

sol =  solvesdp(iff(d1,d2),d1)
assert(sol.problem == 0);
assert(abs(double(d2)-0) <= 1e-3)
