function tests = test_operator_implies
tests = functiontests(localfunctions);

function test1(testCase)

% Binary variable implies LP constriants
sdpvar y u
binvar x
sol = optimize([-10<=u<=10,implies(x,u>=3+2*x)],-x);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(x)-1) <= 1e-4);
testCase.assertTrue(value(u)>4.999);

function test2(testCase)

% Binary variable implies LP constriants
sdpvar y u
binvar x
sol = optimize([-10<=u<=10,-6 <= y <= 6,implies(x,[u>=3+2*x, y >= u])],-x);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(x)-1) <= 1e-4);
testCase.assertTrue(value(u)>4.999);
testCase.assertTrue(value(y-u)>-0.0001);

function test3(testCase)

% Binary variable implies LP constriants
sdpvar y u
binvar x
sol = optimize([-10<=u<=10,-4 <= y <= 4,implies(x,[u>=3+2*x, y >= u])],-x);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(x)-0) <= 1e-4);
testCase.assertTrue(value(u)>=-10);
testCase.assertTrue(value(y)>=-4);

function test4(testCase)
binvar x
sdpvar u y
sol = optimize([-10<=u<=10,-4 <= y <= 4,implies(x,[u==3+2*x])],-x);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(x)-1) <= 1e-4);
testCase.assertTrue(abs(value(u)-5) <= 1e-4);

function test5(testCase)
binvar x
sdpvar u y
sol = optimize([-10<=u<=10,-4 <= y <= 4,implies(x,[u==3+2*x, y >= 2])],-x);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(x)-1) <= 1e-4);
testCase.assertTrue(abs(value(u)-5) <= 1e-4);
testCase.assertTrue(value(y)>=2);

function test6(testCase)
binvar x
sdpvar u y
sol = optimize([-10<=u<=10,-4 <= y <= 4,implies(1-x, u==0.3),implies(x,[u==3+2*x, y >= 6])],-x);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(x)-0) <= 1e-4);
testCase.assertTrue(abs(value(u)-0.3) <= 1e-4);

function test7(testCase)
sdpvar x u
sol = optimize([-10<=u<=10,-1<=x<=1,implies(x>=0,u==-x),implies(x<=0,u==-2*x)],x);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(x)--1) <= 1e-4);
testCase.assertTrue(abs(value(u)-2) <= 1e-4);

function test8(testCase)
sdpvar x u
sol = optimize([-10<=u<=10,-1<=x<=1,implies(x>=0,u==-x),implies(x<=0,u==-2*x)],-x);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(x)-1) <= 1e-4);
testCase.assertTrue(abs(value(u)--1) <= 1e-4);

function test9(testCase)
sdpvar x u
sol = optimize([-.5<=u<=.5,-1<=x<=1,implies(x>=0,u==-x),implies(x<=0,u==-2*x)],-x);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(x)-.5) <= 1e-4);
testCase.assertTrue(abs(value(u)--.5) <= 1e-4);

function test11(testCase)
sdpvar x u
ops = sdpsettings;
x = sdpvar(2,1);
sol = optimize([-10<=u<=10,-1<=x<=1,implies(x>=0,u==-x(1)),implies(x<=0,u==-2*x(1))],(x+.1)'*(x+.1),ops)
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(u)-.2) <= 1e-4);

sol = optimize([-10<=u<=10,-1<=x<=1,implies(x>=0,u==-x(1)),implies(x<=0,u==-2*x(1))],(x-.1)'*(x-.1),ops);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(u)--.1) <= 1e-4);

sol = optimize([-10<=u<=10,-1<=x<=1,implies(x>=0,u==-x(1)),implies(x<=0,u==-2*x(1))],(x-[.1;-.1])'*(x-[.1;-.1])+u^2,ops);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(norm(value(x)-[.1;-0.1]) <= 1e-4);
testCase.assertTrue(abs(value(u)-0) <= 1e-4);

function test12(testCase)
ops = sdpsettings('quadprog.Algorithm','interior-point-convex');
x = sdpvar(2,1);
u = sdpvar(1);
optimize([-10 <= [x;u] <= 10,implies(x==0.5,u == pi)],(x-.5)'*(x-.5),ops)
testCase.assertTrue(norm(value(x)-.5) <= 1e-3);

optimize([-10 <= [x;u] <= 10,implies(x==0.5,u == 15,1e-3)],(x-.5)'*(x-.5),ops)
testCase.assertTrue(~any(abs(value(x)-0.5)>1e-3));


