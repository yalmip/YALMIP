function test_misc_normconvert

sdpvar x(2,1);
sol = solvesdp([x>=0,norm(x)<= 1],sum(x),sdpsettings('solver','fmincon'))
assertTrue(sol.problem == 0);