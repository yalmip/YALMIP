function tests = test_operator_entropy
tests = functiontests(localfunctions);

function test1(testCase)
x = sdpvar(2,1);
assign(x,1);
sol = optimize((x >= 0.1),-entropy(x),sdpsettings('solver','fmincon'));

testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(sum(x))-0.73) <= 1e-2);
testCase.assertTrue(abs(-entropy(value(x))--0.7357588) <= 1e-3);

function test2(testCase)
x1 = sdpvar(1,1);
y = sdpvar(1,1);
x2 = sdpvar(1,1);
x = [x1;x2];
sol = optimize(([x;y] >= 0.1),-entropy(x)+y,sdpsettings('solver','fmincon'));

testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(sum(x))-0.73) <= 1e-2);

