function tests = test_operator_optimizer19
tests = functiontests(localfunctions);

function test1(dummy)
if 0
    yalmip('clear')
    sdpvar x a b c y
    optimize([0 <= [x y] <= 10, x + y <= 5],(x-3)^2+5*(y-4)^2,sdpsettings('verbose',0));
    sol = value([x y])
    assert(norm(sol-[1 + 1/3;3 + 2/3]') <= 1e-2);
end

sdpvar x a b c y

ops = sdpsettings('solver','gurobi');
P = optimizer([0 <= [x y] <= 10, x + y <= c],(x-a)^2+5*(y-b)^2,sdpsettings('solver','gurobi'),{a,b,c},[x y]);
PA = P{{3,[],[]}};
PB = PA{{4,[]}};
sol = PB{{5}}
assert(norm(sol-[1 + 1/3;3 + 2/3]') <= 1e-2);
PC = PB{{5},'nosolve'};
sol = PC{[]}
assert(norm(sol-[1 + 1/3;3 + 2/3]') <= 1e-2);

P = optimizer([0 <= [x y] <= 10, x + y <= c],(x-a)^2+5*(y-b)^2,ops,{c,b,a},[x y]);
PC = P{{5,[],[]}};
PB = PC{{4,[]}};
sol = PB{{3}}
assert(norm(sol-[1 + 1/3;3 + 2/3]') <= 1e-2);
PA = PB{{3},'nosolve'};
sol = PA{[]}
assert(norm(sol-[1 + 1/3;3 + 2/3]') <= 1e-2);

P = optimizer([0 <= [x y] <= 10, x + y <= c],(x-a)^2+5*(y-b)^2,ops,{b,c,a},[x y]);
PC = P{{[4],[],[]}};
PB = PC{{5,[]}};
sol = PB{{3}}
assert(norm(sol-[1 + 1/3;3 + 2/3]') <= 1e-2);
PA = PB{{3},'nosolve'};
sol = PA{[]}
assert(norm(sol-[1 + 1/3;3 + 2/3]') <= 1e-2);

P = optimizer([0 <= [x y] <= 10, x + y <= c],(x-a)^2+5*(y-b)^2,ops,{a,c,b},[x y]);
PC = P{{[3],[],[]}};
PB = PC{{5,[]}};
sol = PB{{4}}
assert(norm(sol-[1 + 1/3;3 + 2/3]') <= 1e-2);
PA = PB{{4},'nosolve'};
sol = PA{[]}
assert(norm(sol-[1 + 1/3;3 + 2/3]') <= 1e-2);

sdpvar x b c y a
P = optimizer([0 <= [x y] <= 10, x + y <= c],(x-a)^2+5*(y-b)^2,ops,{a,c,b},[x y]);
PC = P{{[3],[],[]}};
PB = PC{{5,[]}};
sol = PB{{4}}
assert(norm(sol-[1 + 1/3;3 + 2/3]') <= 1e-2);

PA = PB{{4},'nosolve'};
sol = PA{[]}
assert(norm(sol-[1 + 1/3;3 + 2/3]') <= 1e-2);

sdpvar b x c y a
P = optimizer([0 <= [x y] <= 10, x + y <= c],(x-a)^2+5*(y-b)^2,ops,{a,c,b},[x y]);
PC = P{{[3],[],[]}};
PB = PC{{5,[]}};
PA = PB{{4}}
PA = PB{{4},'nosolve'};
sol = PA{[]}
assert(norm(sol-[1 + 1/3;3 + 2/3]') <= 1e-2);


sdpvar b x c y a
P = optimizer([0 <= [x y] <= 10, x + y <= c],(x-a)^2+5*(y-b)^2,ops,{a,c,b},[x y]);
PC = P{{[3],[],[]}};
PB = PC{{5,[]}};
PA = PB{{4},'nosolve'};
sol = PA{[]}
assert(norm(sol-[1 + 1/3;3 + 2/3]') <= 1e-2);


sdpvar b x c y a
P = optimizer([0 <= [x y] <= 10, x + y <= c],(x-a)^2+5*(y-b)^2,ops,{a,c,b},[x y]);
PC = P{{[3],[5],[]}};
PA = PC{{4},'nosolve'};
sol = PA{[]}
assert(norm(sol-[1 + 1/3;3 + 2/3]') <= 1e-2);


sdpvar b x c y a
P = optimizer([0 <= [x y] <= 10, x + y <= c],(x-a)^2+5*(y-b)^2,ops,{a,c,b},[x y]);
PC = P{{[3],[],[]}};
PA = PC{{5,4},'nosolve'};
sol = PA{[]}
assert(norm(sol-[1 + 1/3;3 + 2/3]') <= 1e-2);


sdpvar b x c y a
P = optimizer([0 <= [x y] <= 10, x + y <= c],(x-a)^2+5*(y-b)^2,ops,{a,c,b},[x y]);
PC = P{{[3],[5],[]}};
PA = PC{{4},'nosolve'};
sol = PA{[]}
assert(norm(sol-[1 + 1/3;3 + 2/3]') <= 1e-2);


P = optimizer([0 <= [x y] <= 10, x + y <= c],(x-a)^2+5*(y-b)^2,ops,{a,c,b},[x y]);
PA = P{[a == 3]};
PB = PA{[b == 4]};
sol = PB{[c == 5]}

P = optimizer([0 <= [x y] <= 10, x + y <= c],(x-a)^2+5*(y-b)^2,ops,{a,c,b},[x y]);
PC = P{[c == 5]};
sol = PC{[a == 3, -b == -4]}
assert(norm(sol-[1 + 1/3;3 + 2/3]') <= 1e-2);

yalmip('clear')
sdpvar x y a b c
ops = sdpsettings('solver','gurobi');
P = optimizer([0 <= [x y] <= 10, x + y <= exp(c)],(x-a)^2+5*(y-exp(b))^2,ops,{a,b,c},[x y]);
PC = P{[a == 3]}
sol = PC{[c == log(5), b == log(4)]}
assert(norm(sol-[1 + 1/3;3 + 2/3]') <= 1e-2);

yalmip('clear')
sdpvar x y a b c
ops = sdpsettings('solver','gurobi');
P = optimizer([0 <= [x y] <= 10, x + y <= exp(c)],(x-a)^2+5*(y-exp(b))^2,ops,{a,b,c},[x y]);
PC = P{[a == 3,c == log(5)]}
sol = PC{[b == log(4)]}
assert(norm(sol-[1 + 1/3;3 + 2/3]') <= 1e-2);