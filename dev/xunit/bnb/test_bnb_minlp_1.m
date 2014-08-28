function test_bnb_minlp_1

randn('seed',12345);
rand('seed',12345);

x = sdpvar(5,1);
A = randn(15,5);
b = rand(15,1)*10;

obj = sum(x) + sum((x-3).^4);
constraints = (A*x <= b) + (integer(x));
sol = solvesdp(constraints,obj,sdpsettings('bnb.solver','fmincon','warning',0))

mbg_asserttrue(sol.problem == 0);
