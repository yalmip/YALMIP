function test_abs

sdpvar x y
obj = abs(1+abs(x-5));
sol = solvesdp([],obj)

mbg_asserttrue(sol.problem == 0)
mbg_asserttolequal(double(obj),1, 1e-4);

sdpvar x y
F = set(abs(1+abs(x-5)) + abs(y)<=10) 
obj = -x
sol = solvesdp(F,obj)
mbg_asserttolequal(double(obj),-14, 1e-4);

