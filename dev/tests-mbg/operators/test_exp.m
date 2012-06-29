function test_exp


sdpvar x y
obj = -x;
sol = solvesdp(set(exp(2*x + 1) <= 3),obj,sdpsettings('solver','fmincon'));

mbg_asserttrue(sol.problem == 0)
mbg_asserttolequal(double(obj),-0.04930614506222, 1e-4);

sdpvar x y
obj = -x-y;
sol = solvesdp(set(exp(max([2*x+1 3*y+2])) <= 3),obj,sdpsettings('solver','fmincon'));

mbg_asserttrue(sol.problem == 0)
mbg_asserttolequal(double(obj),0.25115642610991, 1e-4);


sdpvar x y
obj = -x-y;
sol = solvesdp(set(exp(min([2*x+1 3*y+2])) <= 3),obj,sdpsettings('solver','fmincon','warning',0));
mbg_asserttrue(sol.problem == -4)

sdpvar x y
obj = -x-y;
sol = solvesdp(set(max([exp(2*x+1) exp(3*y+2)]) <= 3),obj,sdpsettings('solver','fmincon'));

mbg_asserttrue(sol.problem == 0)
mbg_asserttolequal(double(obj),0.25115642610991, 1e-4);
