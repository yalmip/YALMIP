function tests = test_misc_normconvert
tests = functiontests(localfunctions);

function test1(dummy)
% Tests bug #282
yalmip('clear')
sdpvar theta
d = cos(theta);
e = sdpvar(1);
f = sdpvar(1);
sol = optimize([e == d,e==e*f],e^2)
assert(sol.problem == 0);