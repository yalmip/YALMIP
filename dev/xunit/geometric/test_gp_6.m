function test_gp_6

q = sdpvar(1,1);
F = (q >= 0);
obj = (1+q)^2.5,
sol = solvesdp(F,obj,sdpsettings('debug',1,'solver','fmincon-geometric'));

mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(double(obj),1,1e-5);