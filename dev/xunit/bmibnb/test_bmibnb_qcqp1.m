function qcqp1

x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);

obj = -2*x1+x2-x3;

F = set(x1*(4*x1-4*x2+4*x3-20)+x2*(2*x2-2*x3+9)+x3*(2*x3-13)+24>=0);
F = F + set(4-(x1+x2+x3)>=0);
F = F + set(6-(3*x2+x3)>=0);
F = F + set(x1>=0);
F = F + set(2-x1>=0);
F = F + set(x2>=0);
F = F + set(x3>=0);
F = F + set(3-x3>=0);

sol = solvesdp(F,obj,sdpsettings('solver','bmibnb'))

mbg_asserttolequal(double(obj),-4, 1e-3);