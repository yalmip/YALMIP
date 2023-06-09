function tests = test_sos_matrix_2
tests = functiontests(localfunctions);

function test1(testCase)

sdpvar x1 x2
P = [1+x1^2 -x1+x2+x1^2;-x1+x2+x1^2 2*x1^2-2*x1*x2+x2^2];
[sol,v,Q] = solvesos((sos(P)),[],sdpsettings('verbose',0,'sedumi.free',0));

diff = clean(v{1}'*Q{1}*v{1}-P,1e-6);

testCase.assertTrue(isequal(diff,[0 0;0 0]))
