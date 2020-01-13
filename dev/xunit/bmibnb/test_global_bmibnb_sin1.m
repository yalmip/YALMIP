function tests = test_global_bmibnb_sin1
tests = functiontests(localfunctions);

function test1(dummy)

sdpvar x
obj = sin(10*x)+abs(x)
sol = optimize((-pi <= x <= pi),obj,sdpsettings('solver','bmibnb'));

assert(sol.problem==0)
assert(abs(value(obj)--0.84792) <= 1e-4)