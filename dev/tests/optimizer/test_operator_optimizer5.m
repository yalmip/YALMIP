function tests = test_operator_optimizer5
tests = functiontests(localfunctions);

function test1(testCase)
sdpvar x y u z 
A = sdpvar(2,3);
P1 = optimizer([x <= u,y <= z], -x-y,[],[u;z],[x;y]);
sol1 = P1{[2;3]};
testCase.assertTrue(norm(sol1-[2;3]) <= 1e-4);
P1 = optimizer([x <= u,y <= z], -x-y,[],[z;u],[x;y]);
sol1 = P1{[3;2]};
testCase.assertTrue(norm(sol1-[2;3]) <= 1e-4);
P1 = optimizer([x <= u,y <= z], -x-y,[],[z;u],[y;x]);
sol1 = P1{[3;2]};
testCase.assertTrue(norm(sol1-[3;2]) <= 1e-4);
P1 = optimizer([x <= u,y <= z], -x-y,[],[u z],[x;y]);
sol1 = P1{[2 3]};
testCase.assertTrue(norm(sol1-[2;3]) <= 1e-4);
P1 = optimizer([x <= u,y <= z], -x-y,[],[u z],[x^2;y]);
sol1 = P1{[2 3]};
testCase.assertTrue(norm(sol1-[4;3]) <= 1e-4);
P1 = optimizer([x <= u,y <= z], -x-y,[],[u z],[x^2+4*y+4 y]);
sol1 = P1{[2 3]};
testCase.assertTrue(norm(sol1-[20 3]) <= 1e-4);
P1 = optimizer([x <= u,y <= z], -x-y,[],[u z u],[x^2+4*y+4;y]);
sol1 = P1{[2 3 2]};
testCase.assertTrue(norm(sol1-[20;3]) <= 1e-4);
P1 = optimizer([x <= u,y <= z], -x-y,[],[z z u]',[x^2+4*y+4;y]);
sol1 = P1{[3 3 2]'};
testCase.assertTrue(norm(sol1-[20;3]) <= 1e-4);
P1 = optimizer([x <= u,y <= z], -x-y,[],[z z;u u]',[x^2+4*y+4;y]);
sol1 = P1{[3 3;2 2]'};
testCase.assertTrue(norm(sol1-[20;3]) <= 1e-4);
P1 = optimizer([x <= u,y <= z], -x-y,[],[u z],[x^2+4*y+4;y]);
sol1 = P1{[2 3]};
testCase.assertTrue(norm(sol1-[20;3]) <= 1e-4);
P3 = optimizer([x <= u,y <= z], -x-y,[],{u,z},{x,y});
sol3 = P3{{4,5}};
testCase.assertTrue(isa(sol3,'cell') && length(sol3)==2 && abs(sol3{2}-5) <= 1e-5);
P4 = optimizer([x <= u,y <= z],norm(A-x)-y,[],{u,z,A},{x,y,A});
sol4 = P4{{4,5,ones(2,3)}};
testCase.assertTrue(isa(sol4,'cell') & length(sol4)==3);