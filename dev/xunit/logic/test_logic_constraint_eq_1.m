function test_logic_constraint_eq_1

binvar a b
intvar x
F = [ a == (b & (x<=4.5)), -15 <= x <= 15]
sol = solvesdp(F,a^2);
mbg_asserttrue(sol.problem == 0);
% Since a is 0, b has to be false or constraint violated
mbg_asserttrue(double(a)==0 && (double(b)==0 || double(x)>=5));

dX = binvar(1);
dY = binvar(1);
sys = [iff(dX,a), iff(dY,[b & (x<=4.5)]), dX == dY];
F = [sys, -15 <= x <= 15]
sol = solvesdp(F,x^2 + 1-b)
mbg_asserttrue(double(a)==1 && (double(b)==1 || double(x)==0));