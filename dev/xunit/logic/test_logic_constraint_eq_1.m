function tests = test_logic_constraint_eq_1
tests = functiontests(localfunctions);

function test1(dummy)
binvar a b
intvar x
F = [ a == (b & (x<=4.5)), -15 <= x <= 15]
sol = optimize(F,a^2);
assert(sol.problem == 0);
% Since a is 0, b has to be false or constraint violated
assert(value(a)==0 && (value(b)==0 || value(x)>=5));

F = true(or(0,a));
sol = optimize(F,a);
assert(sol.problem == 0);
assert(value(a)==1);

F = true(and(0,a));
sol = optimize(F,a);
assert(sol.problem == 1);

F = true(xor(0,a));
sol = optimize(F,a);
assert(sol.problem == 0);
assert(value(a)==1);

F = true(xor(1,a));
sol = optimize(F,-a);
assert(sol.problem == 0);
assert(value(a)==0);

F = [ a == (b & (x<=4.5)), -15 <= x <= 15]
sol = optimize(F,x^2 + 1-b)
assert(value(a)==1 && (value(b)==1 || value(x)==0));

dX = binvar(1);
dY = binvar(1);
sys = [iff(dX,a), iff(dY,[b & (x<=4.5)]), dX == dY];
F = [sys, -15 <= x <= 15]
sol = optimize(F,x^2 + 1-b)
assert(value(a)==1 && (value(b)==1 || value(x)==0));

a = binvar(3,2);
b = binvar(3,2);
optimize([b(1)==0, true(or(a,b))],sum(sum(a-b)))
assert(value(sum(sum(a-b)))==-4);

optimize([b(1)==1, true(and(a,b))],sum(sum(a-b)))
assert(value(sum(sum(a-b)))==0);

optimize([b(1)==1, true(xor(a,b))],sum(sum(-a-b)))
assert(value(sum(sum(-a-b)))==-6);

optimize([true(a & (~b))],sum(sum(a-b)))
assert(value(sum(sum(a-b)))==6);