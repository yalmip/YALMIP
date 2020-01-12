function tests = test_operator_optimizer12
tests = functiontests(localfunctions);

function test1(dummy)

% Tests bug #162 that made eliminatevariables flawed.

sdpvar x y
P = optimizer(cone([1;x])+cone([y;x])+[1 x*y;x*y x+y],x,sdpsettings('solver','+mosek'),y,x)
assert(abs(P{0}) <= 1e-3)
