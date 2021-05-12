function tests = test_diag
tests = functiontests(localfunctions);

function test1(dummy)
P = sdpvar(3,3,'skew');
assert(all([0;0;0] == diag(P)));