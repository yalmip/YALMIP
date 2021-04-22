function tests = test_regress_weird4
tests = functiontests(localfunctions);

function test1(dummy)

% Issue #188
yalmip('clear')
sdpvar x y 
p = sin(1+y*x)^2+cos(y*x);
sdpvar z w
sol = optimize([-1 <= [x y z w] <= 1, p <= 3],p,sdpsettings('solver','bmibnb'));

assert(sol.problem == 0)
assert(abs(value(p)-.5403) <= 1e-2) 