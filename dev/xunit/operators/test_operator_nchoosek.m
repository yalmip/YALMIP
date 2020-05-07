function tests = test_operator_nchoosek
tests = functiontests(localfunctions);

function test1(dummy)
yalmip('clear');
X = sdpvar(10,1);
k = 3;
alpha = rand(10,1).^5;
Model = [sum(X)==40,3<=X<=15];
sol = optimize(Model,sum(alpha.*nchoosek(X,3)));
assert(sol.problem == 0);