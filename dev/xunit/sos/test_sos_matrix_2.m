function test_sos_matrix_2

sdpvar x1 x2
P = [1+x1^2 -x1+x2+x1^2;-x1+x2+x1^2 2*x1^2-2*x1*x2+x2^2];
[sol,v,Q] = solvesos((sos(P)));

diff = clean(v{1}'*Q{1}*v{1}-P,1e-6);

mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(diff,[0 0;0 0]);
