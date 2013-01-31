function test_optimizer4
sdpvar x y u z 
A = sdpvar(2,3);

P1 = optimizer([x <= u,y <= z], -x-y,[],[u;z],[x;y]);

sol1 = P1{[2;3]};

mbg_asserttolequal(sol1,[2;3], 1e-4);

P2 = optimizer([x <= u,y <= z], -x-y,[],{u,z},[x;y]);

%sol2 = P2{[2;3]};

P3 = optimizer([x <= u,y <= z], -x-y,[],{u,z},{x,y});
sol3 = P3{{4,5}};
mbg_asserttrue(isa(sol3,'cell') & length(sol3)==2);

P4 = optimizer([x <= u,y <= z],norm(A-x)-y,[],{u,z,A},{x,y,A});
sol4 = P4{{4,5,ones(2,3)}};
mbg_asserttrue(isa(sol4,'cell') & length(sol4)==3);