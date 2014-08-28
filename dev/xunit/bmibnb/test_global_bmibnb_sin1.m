function sin1

sdpvar x
obj = sin(10*x)+abs(x)
sol = solvesdp((-pi <= x <= pi),obj,sdpsettings('solver','bmibnb'));

mbg_asserttrue(sol.problem==0);
mbg_asserttolequal(double(obj),-0.84792, 1e-4);