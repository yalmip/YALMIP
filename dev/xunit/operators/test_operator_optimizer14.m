function test_operator_optimizer14

sdpvar x a
P = optimizer(x/(1+a^2) >= 1,x^2,sdpsettings('solver','quadprog','verbose',2),a,x);
mbg_asserttolequal(P{7},50,1e-4);