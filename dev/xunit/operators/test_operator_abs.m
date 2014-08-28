function test_operator_abs

sdpvar x y
obj = abs(1+abs(x-5));
sol = solvesdp([],obj)

assertTrue(sol.problem == 0)
assertElementsAlmostEqual(double(obj),1,'absolute', 1e-4);

sdpvar x y
F = (abs(1+abs(x-5)) + abs(y)<=10) 
obj = -x
sol = solvesdp(F,obj)
assertElementsAlmostEqual(double(obj),-14, 'absolute',1e-4);

