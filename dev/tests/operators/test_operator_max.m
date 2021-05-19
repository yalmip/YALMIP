function tests = test_operator_max
tests = functiontests(localfunctions);

function test1(testCase)
% Test for bug #345
a=sdpvar(1,1);
b=sdpvar(1,1);
C=[];
C=[C,a==-1];
C=[C,b==max(a,0)];
sol = optimize(C)
testCase.assertTrue(abs(value(b)) <= 1e-4);

function test2(testCase)
% Test for feature #313
y = (1:3)';
x = intvar(3,1);
[val,loc] = max(x);
sol = optimize([sort(x) == y, loc == 3],x'*x)
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(x(3))-3) <= 1e-4);

[val,loc] = min(x);
sol = optimize([sort(x) == y, loc == 3],x'*x)
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(x(3))-1) <= 1e-4);

function test3(testCase)
sdpvar x y
obj = -x-y;
F = (max([x^2+y^2 x+y]) <= 3);
sol = optimize((max([x^2+y^2 x+y]) <= 3),obj)
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(obj)--sqrt(3/2)*2) <= 1e-4);

obj = -x-y;
sol = optimize((0 <= min([-x^2-y^2 -x-y]) +3),obj)
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(obj)--sqrt(3/2)*2) <= 1e-4);

obj = -x-y;
sol = optimize((max([1 x y x^2]) <= min([-x^2-y^2 -x-y]) +3),obj)
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(obj)--2) <= 1e-4);

obj = -x-y;
sol = optimize((abs([1 x y x^2]) <= min([-x^2-y^2 -x-y]) +3),obj,sdpsettings())
testCase.assertTrue(sol.problem == 0)

obj = -x-y;
sol = optimize((abs([1 x y x^2]) <= min([-x^2-y^2 -x-y]) +3),obj,sdpsettings('allownonconvex',0))
testCase.assertTrue(sol.problem == 14)
