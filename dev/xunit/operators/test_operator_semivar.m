function tests = test_operator_semivar
tests = functiontests(localfunctions);

function test1(dummy)
A = magic(10);
b = A(:,end);
A = A(:,1:5);
x = semivar(5,1);
e = b-A*(x-1);
obj = norm(e,1);
sol = optimize([1 <= x <= 2],obj);

assert(sol.problem == 0)
assert(abs(value(obj)-25.211) <= 1e-3);

% Still not working
obj = e'*e
sol = optimize([1 <= x <= 2],obj);
assert(sol.problem == 0)
assert(abs(value(obj) - 133.2742) <= 1e-3);

