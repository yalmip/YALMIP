function test_lmirank

A = [0.2 0 1 0;0 0.2 0 1;-1 1 0.2 0;1 -1 0 0.2];
B = [0 0 1 0]';
C = [0 1 0 0];
 

X = sdpvar(4,4);
Y = sdpvar(4,4);
 
Bp = null(B')';
Cp = null(C)';
W  = eye(3)*1e-6;
F = set(Bp*(A*X+X*A')*Bp' <= -W) + set(Cp*(Y*A+A'*Y)*Cp' <= -W);
 
F = F + set([X eye(4);eye(4) Y] >= 0);
F = F + set(rank([X eye(4);eye(4) Y]) <= 4+2);
 
sol = solvesdp(F,[],sdpsettings('lmirank.solver','sedumi','sedumi.eps',0))
if sol.problem ~=-2
e = abs(eig(double([X eye(4);eye(4) Y])));
mbg_asserttrue(sol.problem == 0);
mbg_asserttrue(nnz(e > 1e-6) == 6);
end


