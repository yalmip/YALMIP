function tests = test_global_abel
tests = functiontests(localfunctions);

function test1(dummy)

% Issue #187
yalmip('clear')
sdpvar x y
p = sin(1+y)^2+cos(y*x);
sol = optimize([-1 <= [x y] <= 1, p <= 3],p,sdpsettings('solver','bmibnb'));
assert(sol.problem == 0)
assert(abs(value(p)-.5403) <= 1e-2) 