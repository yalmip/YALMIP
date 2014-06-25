function quartic

sdpvar x y
F = set(x^3+y^5 <= 5) + set(y >= 0);
F = F + set(-100 <= [x y] <= 100) % Always bounded domain
obj = -x;

sol = solvesdp(F,obj,sdpsettings('solver','bmibnb'))

mbg_asserttolequal(double(obj),-1.71, 1e-4);