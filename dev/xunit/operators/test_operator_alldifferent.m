function test_operator_alldifferent
n = 4;
x = sdpvar(n,1);
sol = solvesdp((1<=x<=n) + (alldifferent(x)),sum(x))

assertTrue(sol.problem == 0)
assertElementsAlmostEqual(sort(double(x)),(1:n)', 'absolute',1e-4);

x = intvar(1,4);
F = (1 <= x <= 4) + (alldifferent(x))
F = F + (0.5 <= x(4) <= 1.5)
F = F + (3.5 <= x(3) <= 4.5)
sol = solvesdp(F);
assertTrue(sol.problem == 0)
assertElementsAlmostEqual(double(x(4)),1,'absolute', 1e-4);
