function tests = test_sos_matrix_2
tests = functiontests(localfunctions);

function test1(dummy)

sdpvar x1 x2
P = [1+x1^2 -x1+x2+x1^2;-x1+x2+x1^2 2*x1^2-2*x1*x2+x2^2];
[sol,v,Q] = solvesos((sos(P)));

diff = clean(v{1}'*Q{1}*v{1}-P,1e-6);

assert(sol.problem == 0)
assert(isequal(diff,[0 0;0 0]))
