function tests = test_global_bmibnb_sin2
tests = functiontests(localfunctions);

function test1(dummy)

% Tests mixed monomial evals

sdpvar x
obj = sin(cos(x.^2).^2).^2 + 0.1*(x-2).^2;
sol = optimize([-5 <= x <= 5],obj,sdpsettings('allownon',1,'solver','bmibnb','bmibnb.absgaptol',1e-8,'bmibnb.relgaptol',1e-8))

assert(sol.problem == 0)
assert(abs(value(obj)-0.0022650) <= 1e-4)