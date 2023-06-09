function tests = test_operator_optimizer8
tests = functiontests(localfunctions);

function test1(testCase)

% Make sure vectorized multiple-call works
sdpvar x u
P = optimizer([x <= u],-x,[],u,x)
U = P{[7 8 9]};
testCase.assertTrue(norm([7 8 9]-U)<1e-7);