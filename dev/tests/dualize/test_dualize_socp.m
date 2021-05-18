function tests = test_dualize_socp
tests = functiontests(localfunctions);

function test1(testCase)

X = sdpvar(3,3);
x = sdpvar(3,1);
obj = trace(X)+sum(x);
F = (X>=0) + (cone(x(2:end),x(1))) + (trace(X)==x(1)+2*x(2)+3*x(3)+4)+(X(1,3)==8);

sol1  = optimize(F,obj,sdpsettings('verbose',0));
obj1 = value(obj);
p1   = check(F);

sol2 = optimize(F,obj,sdpsettings('dualize',1,'verbose',0));
obj2 = value(obj);
p2   = check(F);

testCase.assertTrue(sol1.problem == 0);


function test2(testCase)

X = sdpvar(3,3);
x = sdpvar(3,1);
obj = trace(X)+sum(x);
F = (X>=0) + (cone(x(2:end),1+x(1))) + (trace(X)==x(1)+2*x(2)+3*x(3)+4)+(X(1,3)==8);

sol1  = optimize(F,obj,sdpsettings('verbose',0));
obj1 = value(obj);
p1   = check(F);

sol2 = optimize(F,obj,sdpsettings('dualize',1,'verbose',0));
obj2 = value(obj);
p2   = check(F);

testCase.assertTrue(sol1.problem == 0);
testCase.assertTrue(sol2.problem == 0);
testCase.assertTrue(abs(obj1 - obj2) <= 1e-4);
testCase.assertTrue(abs(min(p1))<= 1e-4)
testCase.assertTrue(abs(min(p2))<= 1e-4)
testCase.assertTrue(sol2.problem == 0);
testCase.assertTrue(abs(obj1 - obj2) <= 1e-4);
testCase.assertTrue(abs(min(p1))<= 1e-4)
testCase.assertTrue(abs(min(p2))<= 1e-4)


function test3(testCase)

X = sdpvar(3,3);
x = sdpvar(3,1);
obj = trace(X)+sum(x);
F = (X>=0) + (cone(x(2:end),1+x(1))) + (trace(X)==x(1)+2*x(2)+3*x(3)+4)+(X(1,3)==8);

sol1  = optimize(F,obj,sdpsettings('verbose',0));
obj1 = value(obj);
p1   = check(F);

sol2 = optimize(F,obj,sdpsettings('dualize',1,'verbose',0));
obj2 = value(obj);
p2   = check(F);

testCase.assertTrue(sol1.problem == 0);
testCase.assertTrue(sol2.problem == 0);
testCase.assertTrue(abs(obj1 - obj2) <= 1e-4);
testCase.assertTrue(abs(min(p1))<= 1e-4)
testCase.assertTrue(abs(min(p2))<= 1e-4)


function test4(testCase)

X = sdpvar(3,3);
x = sdpvar(3,1);
obj = trace(X)+sum(x);
F = (X>=0) + (cone(1-x(2:end),1+x(1))) + (trace(X)==x(1)+2*x(2)+3*x(3)+4)+(X(1,3)==8);

sol1  = optimize(F,obj,sdpsettings('verbose',0));
obj1 = value(obj);
p1   = check(F);

sol2 = optimize(F,obj,sdpsettings('dualize',1,'verbose',0));
obj2 = value(obj);
p2   = check(F);

testCase.assertTrue(sol1.problem == 0);
testCase.assertTrue(sol2.problem == 0);
testCase.assertTrue(abs(obj1 - obj2) <= 1e-4);
testCase.assertTrue(abs(min(p1))<= 1e-4)
testCase.assertTrue(abs(min(p2))<= 1e-4)