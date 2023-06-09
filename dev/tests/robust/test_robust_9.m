function tests = test_robust_9
tests = functiontests(localfunctions);

function test1(testCase)

yalmip('clear')
sdpvar x w t1 t2;
F = (abs(x) + w <= 1) + (w == t1*(-0.5) + t2*0.5) + ([t1 t2]>=0)+(t1+t2 == 1);
sol = solverobust(F,-x,sdpsettings('verbose',0),[w;t1;t2]);
testCase.assertTrue(abs(value(x) - 1/2) <= 1e-5);
testCase.assertTrue(sol.problem == 0);

function test2(testCase)
sdpvar x w t
F = (x+sum(w) <= 10);
W = (-1/2 <= w <= 1/2);
objective = (x-5)'*(x-5) + x*w;
sol = solverobust(F + W,objective,sdpsettings('verbose',0),[w])
testCase.assertTrue(abs(value(x) - 4.75) <= 1e-2);
