function sahinidis


x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
F = (0 <= x1 <= 6)+(0 <= x2 <= 4) + (x1*x2 <= 4); 

obj = -x1-x2;

sol = solvesdp(F,obj,sdpsettings('solver','bmibnb'))

mbg_asserttolequal(double(obj),-6.66666666, 1e-6);