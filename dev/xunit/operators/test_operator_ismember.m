function tests = test_operator_ismember
tests = functiontests(localfunctions);

function test1(dummy)
x1 = [2;2];
x2 = [2;-2];
x3 = [-2;-2];
x4 = [-2;2];
x = sdpvar(2,1);
P1 = polytope((-0.5 <= x-x1 <= 0.5));
P2 = polytope((-0.5 <= x-x2 <= 0.5));
P3 = polytope((-0.5 <= x-x3 <= 0.5));
P4 = polytope((-0.5 <= x-x4 <= 0.5));
F = ismember(x,[P1 P2 P3 P4]);
sol = optimize(F,-sum(x))

assert(sol.problem == 0)
assert(norm(value(x)-[2.5;2.5]) <= 1e-4);
