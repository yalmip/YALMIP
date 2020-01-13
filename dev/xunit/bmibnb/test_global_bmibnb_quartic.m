function tests = test_global_bmibnb_quartic
tests = functiontests(localfunctions);

function test1(dummy)

sdpvar x y
F = (x^3+y^5 <= 5) + (y >= 0);
F = F + (-100 <= [x y] <= 100) % Always bounded domain
obj = -x;

sol = optimize(F,obj,sdpsettings('solver','bmibnb'))

assert(abs(value(obj)--1.71) <= 1e-4)