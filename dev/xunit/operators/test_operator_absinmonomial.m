function tests = test_operator_absinmonomial
tests = functiontests(localfunctions);

function test1(dummy)
yalmip('clear')
sdpvar x
sol = optimize([-2 <= x <= 1],x*abs(x),sdpsettings('solver','bmibnb'))
assert(sol.problem == 0)
assert(abs(value(x)--2) <= 1e-4);
