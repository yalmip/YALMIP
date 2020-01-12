function tests = test_operator_optimizer15
tests = functiontests(localfunctions);

function test1(dummy)

yalmip('clear')
sdpvar x y z
P = optimizer([x <= y*z],x^2,sdpsettings('solver','+quadprog'),[y;z],x)
assert(abs(P{[-2;5]}--10) <= 1e-4);