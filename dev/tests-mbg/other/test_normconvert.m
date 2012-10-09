function test_normconvert

sdpvar x(2,1);
sol = solvesdp([x>=0,norm(x)<= 1],sum(x),sdpsettings('solver','fmincon'))
mbg_asserttrue(sol.problem == 0);