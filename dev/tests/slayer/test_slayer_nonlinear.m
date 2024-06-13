function tests = test_slayer_nonlinear
tests = functiontests(localfunctions);

function test_disjoint_hiddenconvex(testCase)

sdpvar x y z
P = [-3 <= [x y z] <= 3, 
     [x^2-1 (x+z);(x+z) 1]>=0,
     [y^2-1 (x+y);(x+y) 1]>=0]
ops = sdpsettings('solver','bmibnb','bmibnb.relgaptol',1e-8,'bmibnb.relgaptol',1e-8);
ops.slayer.algorithm = 0;
sol = optimize(P,x^2+y^2+z^2,ops);
testCase.assertTrue(sol.problem == 0)
ops.slayer.algorithm = 1;
sol = optimize(P,x^2+y^2+z^2,ops);
testCase.assertTrue(sol.problem == 0)

ops.slayer.algorithm = 0;
sol = optimize(P,(x-y)^2,ops);
testCase.assertTrue(sol.problem == 0)
ops.slayer.algorithm = 1;
sol = optimize(P,(x-y)^2,ops);
testCase.assertTrue(sol.problem == 0)

ops.slayer.algorithm = 0;
sol = optimize(P,(x-y)^4,ops);
testCase.assertTrue(sol.problem == 0)
ops.slayer.algorithm = 1;
sol = optimize(P,(x-y)^4,ops);
testCase.assertTrue(sol.problem == 0)

P = [-3 <= [x y] <= 3,      
     [y^2-1 (x+y);(x+y) 1]>=0]
ops.slayer.algorithm = 0;
sol = optimize(P,(exp(x)-y)^4,ops);
testCase.assertTrue(sol.problem == 0)
ops.slayer.algorithm = 1;
sol = optimize(P,(exp(x)-y)^4,ops);
testCase.assertTrue(sol.problem == 0)
