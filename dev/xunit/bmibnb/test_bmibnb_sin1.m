function sin1

sdpvar x
obj = sin(10*x)+abs(x)
sol = solvesdp(set(-pi <= x <= pi),obj,sdpsettings('solver','bmibnb'));

mbg_assertfalse(sol.problem);
mbg_asserttolequal(double(obj),-0.84792, 1e-4);