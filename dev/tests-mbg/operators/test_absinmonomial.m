function test_absinmonomial

yalmip('clear')
sdpvar x
sol = solvesdp([-2 <= x <= 1],x*abs(x),sdpsettings('solver','bmibnb'))
mbg_asserttrue(sol.problem == 0)
mbg_asserttolequal(double(x),-2, 1e-4);
