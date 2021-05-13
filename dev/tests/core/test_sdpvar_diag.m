function tests = test_diag
tests = functiontests(localfunctions);

function test1(testCase)
P = sdpvar(3,3,'skew');
testCase.assertTrue(all([0;0;0] == diag(P)));