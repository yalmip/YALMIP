function test_operator_optimizer17

% Test partial instantiation

sdpvar x a y b
ops = sdpsettings('solver','cplex');
P = optimizer([0 <= [x y] <= 10, x + y <= 1],a^3 + (x-a)^2+5*(y-b)^2,ops,{a,b},[x y]);
H = P{{0.5,[]}}
[sol,infeas] = H{0.9};
mbg_asserttolequal(sol,[1/6 5/6],1e-4);
[sol,infeas] = H{0};
mbg_asserttolequal(sol,[.5 0],1e-4);

H = P{{0,[]}}
[sol,infeas] = H{0};
mbg_asserttolequal(sol,[0 0],1e-4);

P = optimizer([0 <= [x y] <= 10, x + y <= 1],a^3 + (x-a-.5)^2+5*(y-b-.2)^2,ops,{a,b},[x y]);
H = P{{0,[]}}
[sol,infeas] = H{0};
mbg_asserttolequal(sol,[.5 0.2],1e-4);

P = optimizer([0 <= [x y] <= 10, x + y <= 1],a^3 + (x-a)^2+5*(y-b)^2,ops,{b,a},[y x]);
H = P{{0.5,[]}}
[sol,infeas] = H{0.9};
mbg_asserttolequal(sol,[ 0.4333    0.5667],1e-3);

sdpvar x a b c y
P = optimizer([0 <= [x y] <= 10, x + y <= c],a^3 + (x-a)^2+5*(y-b)^2,ops,{a,c,b},[y x]);
sol = P{{4,5,6}}
mbg_asserttolequal(sol,[5 0],1e-3);
H = P{{[],5,[]}};
[sol,infeas] = H{{4,6}}
mbg_asserttolequal(sol,[5 0],1e-3);

sol = P{[a == 5, c == 6, b == 5]};
mbg_asserttolequal(sol,[4+1/3 5/3],1e-3);

H = P{[a == 5, b == 5]};
sol = H{[c == 6]}
mbg_asserttolequal(sol,[4+1/3 5/3],1e-3);

sdpvar x a b c y
P = optimizer([0 <= [x y] <= 10, x + y <= c],(x-b*a)^2+5*(y-b)^2,ops,{a,c,b},[y x]);
sol = P{{4,5,6}}
mbg_asserttolequal(sol,[1.8333 3.1667],1e-3);
H = P{{[],5,[]}};
[sol,infeas] = H{{4,6}}
mbg_asserttolequal(sol,[1.8333 3.1667],1e-3);

sdpvar x a b c y
P = optimizer([0 <= [x y] <= 10, x + y <= c],(x-b*a)^2+c*(y-b)^2,ops,{a,c,b},[y x]);
sol = P{{4,5,6}}
mbg_asserttolequal(sol,[1.8333 3.1667],1e-3);
H = P{{[],5,[]}};
[sol,infeas] = H{{4,6}}
mbg_asserttolequal(sol,[1.8333 3.1667],1e-3);



