function tests = test_operator_sort
tests = functiontests(localfunctions);

function test1(dummy)
x = sdpvar(4,1);
z = sdpvar(4,1);

[y,loc] = sort(x);

w = randn(4,1);

sol = optimize((-100 <= x <= 100)+(z == y),norm(x-w,1));

assert(sol.problem == 0);
assert(norm(sort(w)-value(z)) <= 1e-4);

function test2(dummy)
A = ones(20,5);
b = (1:20)';
x = sdpvar(5,1);
e = b-A*x;
F = (mean(x) == median(x)) + (-100 <= x <= 100);
sol = optimize(F,norm(e,1));
assert(sol.problem == 0);
assert(abs(mean(value(x))-median(value(x)))<1e-4);
assert(abs(norm(value(e),1)-100)<1e-4);
