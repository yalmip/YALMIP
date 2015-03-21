function test_misc_normconvert
% Tests bug #282
yalmip('clear')
sdpvar theta
d = cos(theta);
e = sdpvar(1);
f = sdpvar(1);
sol = optimize([e == d,e==e*f],e^2)
assertTrue(sol.problem == 0);