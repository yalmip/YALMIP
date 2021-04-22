function tests = test_cumsum
tests = functiontests(localfunctions);

function test1(dummy)
dx = sdpvar(1,2);
cs1 = [0 cumsum(dx) 1];
cs2 = [0 dx(1) dx(1)+dx(2) 1];
assert(isequal(cs1-cs2,[0 0 0 0]))
