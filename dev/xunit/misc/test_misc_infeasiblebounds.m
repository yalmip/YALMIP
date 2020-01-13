function tests = test_misc_infeasiblebounds
tests = functiontests(localfunctions);

function test1(dummy)
yalmip('clear')
sdpvar a x

F = [ 0 <= x <= 1; a >= 0; a == -0.5 ];

sol = optimize(F, x^3,sdpsettings('solver','bmibnb'))
assert(sol.problem == 1);

sol = optimize(F, x^3,sdpsettings('solver','bnb'))
assert(sol.problem == 1);


