function tests = test_operator_optimizer4
tests = functiontests(localfunctions);

function test1(testCase)
yalmip('clear')
sdpvar x y u z 
P = optimizer([x <= u,y <= z], -x-y,[],[u;z],[x;y]);
sol = optimize([],(P{[u;z]}-[7;2])'*(P{[u;z]}-[7;2]));
testCase.assertTrue(sol.problem == 0)

function test2(testCase)
yalmip('clear')
sdpvar x y u z 
P = optimizer([x <= u,y <= z], -x-y,[],[u;z],[x;y]);
sol = optimize([],(P{1+[u;z]}-[7;2])'*(P{1+[u;z]}-[7;2]));
testCase.assertTrue(sol.problem == 0)