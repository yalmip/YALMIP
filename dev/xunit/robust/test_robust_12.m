function test_robust_12
% Test the duality filter, in particular derivation of c,B,A
yalmip('clear')

x = sdpvar(3,1)
w = sdpvar(2,1)
A = magic(10);
E = A(:,1:2);
A = A(:,2:4);

sol = solvesdp([uncertain(w),-1<=w<=1,sum(w)<=1.5,(A+E*w*ones(1,3))*x <= 10],sum(x),sdpsettings('robust.lplp','duality'))
mbg_asserttolequal(double(sum(x)),-0.16868, 1e-3);
mbg_asserttrue(sol.problem == 0);

sol = solvesdp([uncertain(w),-1<=w<=1,sum(w)<=1.5,(A+E*w*ones(1,3))*x <= 10],sum(x),sdpsettings('robust.lplp',''))
mbg_asserttolequal(double(sum(x)),-0.16868, 1e-3);
mbg_asserttrue(sol.problem == 0);


