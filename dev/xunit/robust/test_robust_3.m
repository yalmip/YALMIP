function tests = test_robust_3
tests = functiontests(localfunctions);

function test1(dummy)
sdpvar x u

sol = solvesdp(((1+u)*x+x>=0.2) + (0.1 <= u <= 0.3) + (uncertain(u)),x)

assert(sol.problem == 0)
assert(abs(double(x) - 9.523809524e-002) <= 1e-5)
