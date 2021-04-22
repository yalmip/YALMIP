function tests = test_misc_normconvert
tests = functiontests(localfunctions);

function test1(dummy)
sdpvar x(2,1);
sol = optimize([x>=0,norm(x)<= 1],sum(x),sdpsettings('solver','fmincon'))
assert(sol.problem == 0);