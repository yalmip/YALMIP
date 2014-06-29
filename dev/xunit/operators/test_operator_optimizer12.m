function test_operator_optimizer12

% Tests bug #162 that made eliminatevariables flawed.

sdpvar x y
P = optimizer(cone([1;x])+cone([y;x])+[1 x*y;x*y x+y],x,sdpsettings('solver','+sdpt3'),y,x)
mbg_asserttolequal(P{0},0,1e-3)
