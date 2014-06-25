function sahinidis


x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
F = set(0 <= x1 <= 6)+set(0 <= x2 <= 4) + set(x1*x2 <= 4); 

obj = -x1-x2;

sol = solvesdp(F,obj,sdpsettings('solver','bmibnb'))

mbg_asserttolequal(double(obj),-6.66666666, 1e-6);