function tests = test_rcone
tests = functiontests(localfunctions);

function test1(testCase)
x = sdpvar(1);
y = sdpvar(1);
z = sdpvar(3,1);

optimize([z>=1, rcone(z,x,y)],x+y);
testCase.assertTrue(abs(value(x*y)-1.5) <= 1e-4);

function test2(testCase)
x = sdpvar(4,1);
sdpvar z y
sol = optimize([(x-1)'*(x-1) <= y*z,y>=0,z>=0], y + z)
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(all(abs(value(x)-1)<1e-4));
testCase.assertTrue(abs(value(y+z))<1e-4);
