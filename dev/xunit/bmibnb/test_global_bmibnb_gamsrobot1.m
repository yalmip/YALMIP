function test_global_bmibnb_gamsrobot1

x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);
x4 = sdpvar(1,1);
t = sdpvar(1,1);
p = x1 - x2 - x3 - x1*x3 + x1*x4 + x2*x3 - x2*x4;
F = (x1+4*x2 <= 8) + (4*x1+x2 <= 12) + (3*x1+4*x2 <= 12) + (2*x3+x4 <= 8) + (x3+2*x4 <= 8) + (x3+x4 <= 5)+([x1;x2;x3;x4]>=0);
F = F+(p<=t);
obj = t;

sol = solvesdp(F,obj,sdpsettings('solver','bmibnb'))

mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(double(obj),-13, 1e-8);