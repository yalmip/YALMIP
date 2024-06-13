function tests = test_operator_optimizer3
tests = functiontests(localfunctions);

function test1(testCase)
yalmip('clear')
sdpvar x
sdpvar z
sdpvar y
P = optimizer([-1 <= y <= 1,y >= z*z, x <= x+x*z],x^2,sdpsettings('solver',''),[y;z],x);
% Should be infeasible 
[~,err] = P{[0;2]};
testCase.assertTrue(err == 1 | err == 12);
P = optimizer([-1 <= y <= 1,[y z;z 1]>=0, x <= x+x*z],x^2,sdpsettings('solver',''),[y;z],x);
[~,err] = P{[0;2]};
testCase.assertTrue(err == 1 | err == 12);

