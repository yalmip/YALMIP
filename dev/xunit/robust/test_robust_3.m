function test_robust_3

sdpvar x u

sol = solvesdp(((1+u)*x+x>=0.2) + (0.1 <= u <= 0.3) + (uncertain(u)),x)

mbg_asserttolequal(sol.problem, 0, 1e-5);
mbg_asserttolequal(double(x),9.523809524e-002, 1e-5);
