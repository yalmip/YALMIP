function tests = test_bnb_minlp_1
tests = functiontests(localfunctions);

function test1(dummy)
randn('seed',12345);
rand('seed',12345);

x = sdpvar(5,1);
A = randn(15,5);
b = rand(15,1)*10;

obj = sum(x) + sum((x-3).^4);
constraints = (A*x <= b) + (integer(x));
sol = optimize(constraints,obj,sdpsettings('solver','bnb','bnb.solver','fmincon','warning',0))

assert(sol.problem == 0);
