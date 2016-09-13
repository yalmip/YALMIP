function test_operator_rcone_automatic

x = sdpvar(4,1);
sdpvar z y

sol = optimize([(x-1)'*(x-1) <= y*z,y>=0,z>=0], y + z)
assertTrue(sol.problem == 0);
assertTrue(all(abs(value(x)-1)<1e-4));
assertTrue(abs(value(y+z))<1e-4);
