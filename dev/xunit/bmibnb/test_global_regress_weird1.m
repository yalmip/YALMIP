function regress_weird1

sdpvar x
sol = solvesdp((-pi <= x <= pi),2^sin(x),sdpsettings('solver','bmibnb'))
mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(double(x), -pi/2, 1e-3);