function test_optimizer8

% Make sure vectorized multiple-call works
sdpvar x u
P = optimizer([x <= u],-x,[],u,x)
U = P{[7 8 9]};
mbg_asserttrue(isequal([7 8 9],U));