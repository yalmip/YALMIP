function test_operator_nchoosek

yalmip('clear');
X = sdpvar(10,1);
k = 3;
alpha = rand(10,1).^5;
Model = [sum(X)==40,3<=X<=15];
sol = optimize(Model,sum(alpha.*nchoosek(X,3)));
assertTrue(sol.problem == 0);