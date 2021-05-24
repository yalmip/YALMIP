function tests = test_bnb_nlp
tests = functiontests(localfunctions);

function test_minlp(testCase)
randn('seed',12345);
rand('seed',12345);

x = sdpvar(5,1);
A = randn(15,5);
b = rand(15,1)*10;

obj = sum(x) + sum((x-3).^4);
constraints = (A*x <= b) + (integer(x));
sol = optimize(constraints,obj,sdpsettings('solver','bnb','bnb.solver','fmincon','verbose',0));

testCase.assertTrue(sol.problem == 0);
