function tests = test_chance_1
tests = functiontests(localfunctions);

function test1(testCase)
yalmip('clear')
a=sdpvar(1,1);
sdpvar t

Model = [probability(a >= t) >= 0.5,uncertain(a,'normal',0,1)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)) <= 1e-3)

function test2(testCase)
yalmip('clear')
a=sdpvar(1,1);
sdpvar t

Model = [probability(a >= t) >= 0.5,uncertain(a,'normalf',0,1)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)) <= 1e-3)

function test3(testCase)
yalmip('clear')
a=sdpvar(1,1);
sdpvar t

Model = [probability(a >= t) >= 0.95,uncertain(a,'normal',0,4)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)--3.2897) <= 1e-3)

function test4(testCase)
yalmip('clear')
a=sdpvar(1,1);
sdpvar t

Model = [probability(a >= t) >= 0.95,uncertain(a,'normalf',0,2)];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)--3.2897) <= 1e-3)

function test5(testCase)
yalmip('clear')
a=sdpvar(1,1);
sdpvar t

sdpvar R
Model = [probability(a >= t) >= 0.95,uncertain(a,'normalf',0,1+R)];
Model = [Model, uncertain(R),0<=R<=1];
optimize(Model,-t)
testCase.assertTrue(abs(value(t)--3.2897) <= 1e-3)

