function tests = test_expand
tests = functiontests(localfunctions);

function test1(testCase)
% Ensure abs only is modelled once
N = 1;
w = sdpvar(N, 1);
previousPos = rand(N, 1);
constraint = [0.01 <= abs(w), abs(w) == 1, abs(w) == 1];
[~,~,~,m] = export(constraint,[],sdpsettings('solver',''));
testCase.assertTrue(length( m.binary_variables) == 1);

