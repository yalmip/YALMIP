function tests = test_operator_alldifferent
tests = functiontests(localfunctions);

function test1(dummy)
n = 4;
x = sdpvar(n,1);
sol = optimize((1<=x<=n) + (alldifferent(x)),sum(x))

assert(sol.problem == 0)
assert(norm(sort(value(x))-(1:n)') <= 1e-4);

function test2(dummy)
x = intvar(1,4);
F = (1 <= x <= 4) + (alldifferent(x))
F = F + (0.5 <= x(4) <= 1.5)
F = F + (3.5 <= x(3) <= 4.5)
sol = optimize(F);
assert(sol.problem == 0)
assert(abs(value(x(4))-1) <= 1e-4);
