function tests = test_operator_optimizer21
tests = functiontests(localfunctions);

function test1(testCase)

x = sdpvar(2,1);
Objective = -sum(x);
Box = [0 <= x <= 5];
v = sdpvar(3,1);
ops = sdpsettings('solver','');
Solver = optimizer(Box,-x(1)-x(2),ops,[],x);
Cuts = optimizer(v(1:2)'*x <= v(3),[],ops,v,x);
[x,d] = Solver();
testCase.assertTrue(d == 0);
testCase.assertTrue(norm(x-[5;5])<=1e-3);
Solver2 = [Solver, Cuts([2;3;4],'nosolve')]
[x,d] = Solver2();
testCase.assertTrue(d == 0);
testCase.assertTrue(norm(x-[2;0])<=1e-3);
Solver3 = [Solver2, Cuts([2;3;3.5],'nosolve')]
[x,d] = Solver3();
testCase.assertTrue(d == 0);
testCase.assertTrue(norm(x-[1.75;0])<=1e-3);
Solver4 = [Cuts([2;3;3],'nosolve'),Solver3];
[x,d] = Solver4();
testCase.assertTrue(d == 0);
testCase.assertTrue(norm(x-[1.5;0])<=1e-3);


