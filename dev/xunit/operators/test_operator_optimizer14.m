function tests = test_operator_optimizer14
tests = functiontests(localfunctions);

function test1(dummy)

sdpvar x a
P = optimizer(x/(1+a^2) >= 1,x^2,sdpsettings('solver','quadprog','verbose',2),a,x);
assert(abs(P{7}-50) <= 1e-4);