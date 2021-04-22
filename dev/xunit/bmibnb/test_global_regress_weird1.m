function tests = test_global_regress_weird1
tests = functiontests(localfunctions);

function test1(dummy)

sdpvar x
sol = optimize((-pi <= x <= pi),2^sin(x),sdpsettings('solver','bmibnb'))
assert(sol.problem == 0)
assert(abs(value(x)--pi/2) <= 1e-3)
