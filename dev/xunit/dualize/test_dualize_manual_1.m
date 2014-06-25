function test_dualize_manual_1

X = sdpvar(2,2);
t = sdpvar(2,1);
Y = sdpvar(3,3);
obj = trace(X)+trace(Y)+5*sum(t);

F = set(sum(X) == 6+pi*t(1)) + set(diag(Y) == -2+exp(1)*t(2))
F = F + set(Y>=0) + set(X>=0);

sol = solvesdp(F,obj,sdpsettings('dualize',1));

mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(double(t),[-1.90985931710276; 0.73575888234282], 1e-5);
