function gamsgpkl

x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
t = sdpvar(1,1);

F = set([x1;x2]>=0);
F = F + set(- x1 + x2 <= 1);
F = F + set(x1 - x2 <= 1);
F = F + set(- x1 + 2*x2 <= 3);
F = F + set(2*x1 - x2 <= 3);
obj = (2*x1 - 2*x1*x1 + 2*x1*x2 + 3*x2 - 2*x2*x2);
F = F + set(obj<=t);
sol = solvesdp(F,obj,sdpsettings('solver','bmibnb'))

mbg_asserttolequal(double(obj),-3, 1e-4);