function tests = test_logic_mixed
tests = functiontests(localfunctions);

function test1(dummy)
binvar d1 d2

sol =  optimize(~(d1==d2),d2)
assert(sol.problem == 0);
assert(abs(value(d1 + d2)-1) <= 1e-3)

sol =  optimize((d1~=d2),d1+d2)
assert(sol.problem == 0);
assert(abs(value(d1 + d2)-1) <= 1e-3)

sol =  optimize(iff(d1,d2),d1+d2)
assert(sol.problem == 0);
assert(abs(value(d1 + d2)-0) <= 1e-3)

sol =  optimize(iff(d1,d2),d1)
assert(sol.problem == 0);
assert(abs(value(d2)-0) <= 1e-3)
