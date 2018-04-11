function test_gp_2

t1 = sdpvar(1,1);
t2 = sdpvar(1,1);
t3 = sdpvar(1,1);
t = [t1 t2 t3];
obj = (40*t1^-1*t2^-0.5*t3^-1)+(20*t1*t3)+(40*t1*t2*t3);
F = ((1/3)*t1^-2*t2^-2+(4/3)*t2^0.5*t3^-1 <= 1);
F = [F, t >= 0];
sol = solvesdp(F + (t1-t2 <= 1),obj,sdpsettings('solver','fmincon-geometric'));

% mbg_asserttolequal(sol.problem,-8);

