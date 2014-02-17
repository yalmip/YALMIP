function test_cutsdp
c=[1 -2 3]';
x = binvar(3,1);
sdpvar t
A = [-1 2 0;-3 -4 1;0 0 -2];
P = sdpvar(3,3);
a = sdpvar(3,1);
F = [P >= 0, P*A'+A*P <= 0];
F = [F, a<=9*x];
F = [F, a>=-9*x];
F = [F, trace(P) == 1];
F = [F, [t (c+a)';c+a P] >=0];
ops = sdpsettings('solver','cutsdp','cutsdp.maxiter',10);
sol = solvesdp(F,t,ops)
% Just assert that it ran
mbg_asserttrue(sol.problem ~= -1);
