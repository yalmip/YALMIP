function tests = test_sdpvar_replace
tests = functiontests(localfunctions);

function test1(testCase)
sdpvar t
p = 1+t^2;
p2 = replace(p,t,2*t);
testCase.assertTrue(isequal(p2-(1+4*t^2),0))

function test2(testCase)
% Checks that the 0^0 bug in MATLAB6.5 LINUX
% is avoided
sdpvar x t
p = x^2+t;
y = replace(p,t,0);
testCase.assertTrue(isequal(getbase(y), getbase(x^2)))