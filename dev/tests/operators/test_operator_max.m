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
