function test_logic_constraint_eq_1

binvar a b
intvar x
F = [ a == (b & (x<=4.5)), -15 <= x <= 15]
sol = solvesdp(F,a^2);
mbg_asserttrue(sol.problem == 0);
% Since a is 0, b has to be false or constraint violated
mbg_asserttrue(double(a)==0 && (double(b)==0 || double(x)>=5));

F = true(or(0,a));
sol = solvesdp(F,a);
mbg_asserttrue(sol.problem == 0);
mbg_asserttrue(double(a)==1);

F = true(and(0,a));
sol = solvesdp(F,a);
mbg_asserttrue(sol.problem == 1);

F = true(xor(0,a));
sol = solvesdp(F,a);
mbg_asserttrue(sol.problem == 0);
mbg_asserttrue(double(a)==1);

F = true(xor(1,a));
sol = solvesdp(F,-a);
mbg_asserttrue(sol.problem == 0);
mbg_asserttrue(double(a)==0);

F = [ a == (b & (x<=4.5)), -15 <= x <= 15]
sol = solvesdp(F,x^2 + 1-b)
mbg_asserttrue(double(a)==1 && (double(b)==1 || double(x)==0));

dX = binvar(1);
dY = binvar(1);
sys = [iff(dX,a), iff(dY,[b & (x<=4.5)]), dX == dY];
F = [sys, -15 <= x <= 15]
sol = solvesdp(F,x^2 + 1-b)
mbg_asserttrue(double(a)==1 && (double(b)==1 || double(x)==0));

a = binvar(3,2);
b = binvar(3,2);
optimize([b(1)==0, true(or(a,b))],sum(sum(a-b)))
mbg_asserttrue(double(sum(sum(a-b)))==-4);

optimize([b(1)==1, true(and(a,b))],sum(sum(a-b)))
mbg_asserttrue(double(sum(sum(a-b)))==0);

optimize([b(1)==1, true(xor(a,b))],sum(sum(-a-b)))
mbg_asserttrue(double(sum(sum(-a-b)))==-6);

optimize([true(a & (~b))],sum(sum(a-b)))
mbg_asserttrue(double(sum(sum(a-b)))==6);