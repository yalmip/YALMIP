function tests = test_operator_optimizer18
tests = functiontests(localfunctions);

function test1(dummy)
ops = sdpsettings('solver','gurobi');
x = sdpvar(2,1);
c = sdpvar(2,1);
P = optimizer([c'*x <= 1],-sum(x),ops,{c},[x]);
S = [];
for i = 1:500
    ci = randn(2,1);ci = ci/norm(ci);
    S = [S,P{c == ci,'nosolve'}];
end
sol = S{[]};
assert(norm(sol - [.7;.7]) <= 2e-1);
