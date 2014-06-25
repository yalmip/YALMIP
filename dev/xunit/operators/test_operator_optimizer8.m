function test_operator_optimizer8

% Make sure vectorized multiple-call works
sdpvar x u
P = optimizer([x <= u],-x,[],u,x)
U = P{[7 8 9]};
mbg_asserttrue(norm([7 8 9]-U)<1e-7);