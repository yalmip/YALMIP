function test_optimizer4
yalmip('clear')
sdpvar x y u z 
P = optimizer([x <= u,y <= z], -x-y,[],[u;z],[x;y]);
sol = solvesdp([],(P{[u;z]}-[7;2])'*(P{[u;z]}-[7;2]));
mbg_asserttrue(sol.problem == 0)