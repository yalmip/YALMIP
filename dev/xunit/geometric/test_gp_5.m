function test_gp_5

t1 = sdpvar(1,1);
t2 = sdpvar(1,1);
t3 = sdpvar(1,1);
obj = (40*t1^-1*t2^-0.5*t3^-1)+(20*t1*t3)+(40*t1*t2*t3);

F = (max((1/3)*t1^-2*t2^-2+(4/3)*t2^0.5*t3^-1,0.25*t1*t2) <= min((t1+0.5*t2)^-1,t2));
F = F + ((2*t1+3*t2^-1)^0.5 <= 2);
F = [F, [t1 t2 t3] >= 0];

sol = solvesdp(F,obj);

mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(double(obj),2.359439050512407e+002,1e-3);
mbg_asserttolequal(double([t1 t2 t3]), [0.76467168678701   1.23304260692267   4.24155022707061], 1e-3);