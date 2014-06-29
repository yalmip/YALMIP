function test_operator_absinmonomial

yalmip('clear')
sdpvar x
sol = solvesdp([-2 <= x <= 1],x*abs(x),sdpsettings('solver','bmibnb'))
assertTrue(sol.problem == 0)
assertElementsAlmostEqual(double(x),-2, 'absolute',1e-4);
