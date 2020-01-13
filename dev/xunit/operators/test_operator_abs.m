function tests = test_operator_abs
tests = functiontests(localfunctions);

function test1(dummy)

sdpvar x y
obj = abs(1+abs(x-5));
sol = optimize([],obj)

assert(sol.problem == 0)
assert(abs(value(obj) - 1) <= 1e-4);

sdpvar x y
F = (abs(1+abs(x-5)) + abs(y)<=10) 
obj = -x
sol = optimize(F,obj)
assert(abs(value(obj)--14) <= 1e-4);

