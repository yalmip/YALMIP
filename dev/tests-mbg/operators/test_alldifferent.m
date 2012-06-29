function test_alldifferent
n = 4;
x = sdpvar(n,1);
sol = solvesdp(set(1<=x<=n) + set(alldifferent(x)),sum(x))

mbg_asserttrue(sol.problem == 0)
mbg_asserttolequal(sort(double(x)),(1:n)', 1e-4);

x = intvar(1,4);
F = set(1 <= x <= 4) + set(alldifferent(x))
F = F + set(0.5 <= x(4) <= 1.5)
F = F + set(3.5 <= x(3) <= 4.5)
sol = solvesdp(F);
mbg_asserttrue(sol.problem == 0)
mbg_asserttolequal(double(x(4)),1, 1e-4);
