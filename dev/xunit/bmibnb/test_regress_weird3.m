function tests = test_regress_weird3
tests = functiontests(localfunctions);

function test1(dummy)

sdpvar x
obj = sin(sin(x.^2) + x.^2)+0.01*x.^2-sin(x);
sol = optimize([-2*pi <= x <= 2*pi],obj,sdpsettings('allownonconvex',1,'solver','bmibnb','quadprog.Algorithm','interior-point-convex'));

assert(sol.problem == 0)
assert(abs(value(obj)--1.734) <= 1e-3)