function tests = test_robust_5
tests = functiontests(localfunctions);

function test1(dummy)
% Test the duality filter, in particular derivation of c,B,A
yalmip('clear')

x = sdpvar(3,1)
w = sdpvar(2,1)
A = magic(10);
E = A(:,1:2);
A = A(:,2:4);

sol = optimize([uncertain(w),-1<=w<=1,sum(w)<=1.5,(A+E*w*ones(1,3))*x <= 10],sum(x),sdpsettings('robust.lplp','duality'))
assert(abs(value(sum(x))--0.16868) <=  1e-3);
assert(sol.problem == 0);

sol = optimize([uncertain(w),-1<=w<=1,sum(w)<=1.5,(A+E*w*ones(1,3))*x <= 10],sum(x),sdpsettings('robust.lplp',''))
assert(abs(value(sum(x))--0.16868) <=  1e-3);
assert(sol.problem == 0);


