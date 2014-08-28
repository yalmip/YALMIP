function test_operator_exp


sdpvar x y
obj = -x;
sol = solvesdp((exp(2*x + 1) <= 3),obj,sdpsettings('solver','fmincon'));

assertTrue(sol.problem == 0)
assertElementsAlmostEqual(double(obj),-0.04930614506222, 'absolute',1e-4);

sdpvar x y
obj = -x-y;
sol = solvesdp((exp(max([2*x+1 3*y+2])) <= 3),obj,sdpsettings('solver','fmincon'));

assertTrue(sol.problem == 0)
assertElementsAlmostEqual(double(obj),0.25115642610991,'absolute', 1e-4);


sdpvar x y
obj = -x-y;
sol = solvesdp((exp(min([2*x+1 3*y+2])) <= 3),obj,sdpsettings('solver','fmincon','warning',0));
assertTrue(sol.problem == -4)

sdpvar x y
obj = -x-y;
sol = solvesdp((max([exp(2*x+1) exp(3*y+2)]) <= 3),obj,sdpsettings('solver','fmincon'));

assertTrue(sol.problem == 0)
assertElementsAlmostEqual(double(obj),0.25115642610991, 'absolute',1e-4);
