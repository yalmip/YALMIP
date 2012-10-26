function test_infeasiblebounds

yalmip('clear')
sdpvar a x

F = [ 0 <= x <= 1; a >= 0; a == -0.5 ];

sol = solvesdp(F, x^3,sdpsettings('solver','bmibnb'))
mbg_asserttrue(sol.problem == 1);

sol = solvesdp(F, x^3,sdpsettings('solver','bnb'))
mbg_asserttrue(sol.problem == 1);


