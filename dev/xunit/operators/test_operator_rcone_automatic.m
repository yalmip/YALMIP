function tests = test_operator_rcone_automatic
tests = functiontests(localfunctions);

function test1(dummy)
x = sdpvar(4,1);
sdpvar z y

sol = optimize([(x-1)'*(x-1) <= y*z,y>=0,z>=0], y + z)
assert(sol.problem == 0);
assert(all(abs(value(x)-1)<1e-4));
assert(abs(value(y+z))<1e-4);
