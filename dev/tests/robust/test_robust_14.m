function tests = test_robust_14
tests = functiontests(localfunctions);

function test1(testCase)
sdpvar t x 
w = sdpvar(2,1);
ops=sdpsettings('verbose',0,'sedumi.free',0);
M = momentmodel([w'*w == 1, w>=0.5]);
P = sosmodel(sos(1 + w(1)*x + w(2)*x^2 - t),[],ops,[t;w]);
sol = optimize([M,P,uncertain(M)],-t,ops);
testCase.assertTrue(sol.problem == 0 || sol.problem == 4);
testCase.assertTrue(abs(value(t)-.625) <= 1e-3);

