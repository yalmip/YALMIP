function test_operator_optimizer15

yalmip('clear')
sdpvar x y z
P = optimizer([x <= y*z],x^2,sdpsettings('solver','+gurobi'),[y;z],x)
mbg_asserttolequal(P{[-2;5]},-10,1e-4);