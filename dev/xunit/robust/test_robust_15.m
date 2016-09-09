function test_robust_15

x = sdpvar(2,1);
w = sdpvar(2,1);
P = eye(2);
S = [120 140;140 180]/4;

w = sdpvar(2,1);
x = sdpvar(2,1);

sol = optimize([(x+w)'*P*(x+w)<=1,w'*S*w <= 1,uncertain(w)],x(1))

mbg_asserttrue(sol.problem == 0);
mbg_asserttolequal(value(x(1)),-0.4, 1e-3);

M = cone([1;chol(P)*(x+w)])
W = cone([1;chol(S)*w]);
sol = optimize([M,W,uncertain(w)],x(1))

mbg_asserttrue(sol.problem == 0);
mbg_asserttolequal(value(x(1)),-0.4, 1e-3);


