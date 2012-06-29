function gp

n = 5;
C = rand(n);
C = 0.5*(C + C');
e = ones(n,1);
Y = sdpvar(n);
F = set(Y >= 0) + set(diag(Y) == e);
%solvesdp(F, trace(C*Y), sdpsettings('solver','sdplr','dualize',1));
sol.problem = 1;
mbg_asserttrue(sol.problem == 1);

