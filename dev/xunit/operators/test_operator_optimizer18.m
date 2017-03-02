function test_operator_optimizer_18

ops = sdpsettings('solver','cplex');
x = sdpvar(2,1);
c = sdpvar(2,1);
P = optimizer([c'*x <= 1],-sum(x),ops,{c},[x]);
S = [];
for i = 1:500
    ci = randn(2,1);ci = ci/norm(ci);
    S = [S,P{c == ci,'nosolve'}];
end
sol = S{[]};
mbg_asserttolequal(sol,[.7;.7],2e-1);
