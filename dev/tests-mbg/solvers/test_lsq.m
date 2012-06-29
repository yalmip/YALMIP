function test_lsq

n = 5;
m = 10;
C = rand(n,m);
d = rand(n,1);
x = sdpvar(m,1);

sol = solvesdp(x>=0,norm(C*x-d),sdpsettings('solver','lsqnonneg'))
mbg_asserttrue(sol.problem == 0);

sol = solvesdp(x>=0,norm(C*x-d)^2,sdpsettings('solver','lsqnonneg'))
mbg_asserttrue(sol.problem == -4);

sol = solvesdp([sum(x)==1,x>=0],norm(C*x-d),sdpsettings('solver','lsqlin'))
mbg_asserttrue(sol.problem == 0);

x = sdpvar(n,1);
d = rand(m,1);
sol = solvesdp(x>=0,norm(C'*x-d),sdpsettings('solver','lsqnonneg'))
mbg_asserttrue(sol.problem == 0);

sol = solvesdp([sum(x)==1,x>=0],norm(C'*x-d),sdpsettings('solver','lsqlin'))
mbg_asserttrue(sol.problem == 0);


